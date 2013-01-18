#################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

# Author: David Alexander

from __future__ import absolute_import

import math, logging, numpy as np, random
from collections import defaultdict, namedtuple, OrderedDict
from bisect import bisect_left, bisect_right
from itertools import izip
from ..utils import probability_to_qv, noEvidenceConsensusCall
from .. import reference
from ..options import options
from ..Worker import WorkerProcess, WorkerThread
from ..ResultCollector import ResultCollectorProcess, ResultCollectorThread
from ..alignment import AlignmentColumn
from .consumers import consumers, broadcast

#
# The type used in the wire protocol.
#
PluralityLocusSummary = namedtuple("PluralityLocusSummary",
                                   ("referenceId",
                                    "referencePos",
                                    "coverage",
                                    "consensus",
                                    "consensusConfidence",
                                    "consensusFrequency"))

class PluralityAlignmentColumn(AlignmentColumn):
    def __init__(self, *args):
        super(PluralityAlignmentColumn,self).__init__(*args)

    def confidence(self):
        """
        Implements AlignmentColumn.confidence

        Our plurality confidence metric is a phred-transformed
        probability based on the assumption that each read is an
        independent observation of the true base, with error
        probability 0.15.

        For a true genomic sequence G, we have a vector of reads X =
        (X_1, X_2, ..., X_d) with values in S = {s_1, s_2, ..., s_k},
        and corresponding frequencies N = {n_1, n_2, ..., n_k}---where
        we have indexed the sets S and N in descending frequency
        order.

        If we allow that

           Pr(X | G \not\in S) = (0.15)^d;
           Pr(X | G = s_1)     = (0.85)^n_1 * (0.15)^(d-n_1);
           Pr(X | G = s_2)     = (0.85)^n_2 * (0.15)^(d-n_2); etc.

        and then impose an uninformative prior distribution on G, such
        that the events G = s_i \forall i, G \not\in S, are all
        equiprobable (which may or may not be reasonable---revisit
        this later), then the respective posterior probabilities for G
        can be found from the above formulas after normalizing by the
        sum.

        Computations done in the log domain to prevent underflow.
        """
        p, q = 0.85, 0.15
        log_p, log_q = map(math.log, (p, q))
        d = self.coverage()
        unnormalized_log_probs = []
        for si, ni in self.orderedSnippetsAndFrequencies:
            unnormalized_log_probs.append(ni * log_p + (d - ni) * log_q)
        unnormalized_log_probs.append(d * log_q)
        unnormalized_log_probs = np.array(unnormalized_log_probs)

        # calculate the sum in the log domain
        log_sum = np.logaddexp.reduce(unnormalized_log_probs)
        log_probs = unnormalized_log_probs - log_sum
        return probability_to_qv(np.exp(log_probs[0]))

    def consensus(self):
        """Implements AlignmentColumn.consensus"""
        consensusSnippet, consensusFrequency = self.mostFrequentItem()
        return PluralityLocusSummary(self.referenceId,
                                     self.referencePos,
                                     self.coverage(),
                                     consensusSnippet,
                                     self.confidence(),
                                     consensusFrequency)

class PluralityCaller(object):

    # The type of the result that will be passed on to the result
    # collector is:
    #     (ReferenceWindow, list of (Locus, Tuple)]
    # where each Tuple must be convertible to a numpy array with dtype
    # as below:

    def onStart(self):
        random.seed(42)

    def onChunk(self, referenceWindow, alnHits):
        if options.coverage != None:
            coverage = min(len(alnHits), options.coverage)
            alnHits = sorted(random.sample(alnHits, coverage))
        return (referenceWindow, self.plurality(referenceWindow, alnHits))

    @staticmethod
    def plurality(referenceWindow, alignments):
        """
        Convenience method for tabulating plurality in a window.

        Input: referenceWindow, iterable of alignmentHits
        Output: list of (Locus, PluralityLocusSummary)

        Plurality is tabulating the frequency of 'snippets'---strings
        of one or more bases---in the reads, corresponding to each
        reference position.  As the plurality 'consensus', we return
        the most frequent snippet by column.
        """
        refId, refWindowStart, refWindowEnd = referenceWindow
        alignmentColumns = {}

        for hit in alignments:
            alignedRefPositions = hit.referencePositions(orientation="genomic")

            begin = bisect_left(alignedRefPositions, refWindowStart)
            end = bisect_right(alignedRefPositions, refWindowEnd)

            alignedRefPositions = alignedRefPositions[begin:end]
            alignedRef          = hit.reference(orientation="genomic")[begin:end]
            alignedRead         = hit.read(orientation="genomic")[begin:end]

            readBases = []
            for i, (refPos, refBase, readBase) in enumerate(izip(alignedRefPositions,
                                                                 alignedRef,
                                                                 alignedRead)):
                if readBase != "-":
                    readBases.append(readBase)

                if refBase != "-":
                    if (refId, refPos) not in alignmentColumns:
                        alignmentColumns[refId, refPos] = \
                            PluralityAlignmentColumn(refId, refPos, refBase)
                    alignmentColumns[refId, refPos].addReadSnippet("".join(readBases))
                    readBases = []

        # Return (Locus, Tuple)
        return [(locus, alnCol.consensus())
                for locus, alnCol in alignmentColumns.iteritems()]


# define both process and thread-based plurality callers
class PluralityWorkerProcess(PluralityCaller,WorkerProcess):
    pass

class PluralityWorkerThread(PluralityCaller,WorkerThread):
    pass


class PluralityResult(object):
    #
    # The type used in the table built in the result collector process.
    #
    LOCUS_SUMMARY_DTYPE = [ ("referencePos",         np.uint32),
                            ("coverage"            , np.uint32),
                            ("consensus"           , "S8"),
                            ("consensusConfidence" , np.uint8),
                            ("consensusFrequency"  , np.uint32)]

    """
    Collect the plurality results, collate, and output.
    """
    def onStart(self):
        self.consensusByRefId = OrderedDict()
        self.chunksReceivedById = defaultdict(int)
        self.consumers = consumers(options.outputFilenames,
                                   PluralityLocusSummary._fields,
                                   options.variantConfidenceThreshold)

    def onResult(self, result):
        # This is trickier than it really ought to be.  We want to
        # minimize memory consumption, so as soon as we have received
        # all the results for a given reference group, we send the
        # results to disk and then delete its consensus table---which
        # takes up a ton of memory.  In 1.4 we should rewrite things
        # so there is just a top level synchronous loop over refId,
        # and we don't have to worry about this trickiness.
        referenceWindow, pluralityResults = result
        refId = referenceWindow[0]
        refEntry = reference.byId[refId]
        self.chunksReceivedById[refId] += 1

        # Is this a refId we haven't seen before? If so, initialize its
        # consensus table.
        if refId not in self.consensusByRefId:
            self.initTable(refId, refEntry)

        for locus, locusSummary in pluralityResults:
            _, refPos = locus
            assert _ == refId
            self.installInTable(locusSummary, self.consensusByRefId[refId], refPos)

        # Did we finish processing this reference group?  If so, flush to
        # disk!
        if self.chunksReceivedById[refId] == reference.numChunks(refId,
                                                                 options.referenceChunkSize,
                                                                 options.referenceWindow):
            broadcast((refId, self.consensusByRefId[refId]), self.consumers)
            del self.consensusByRefId[refId]

    def onFinish(self):
        logging.info("Analysis completed.")
        for consumer in self.consumers:
            consumer.close()

    def initTable(self, refId, refEntry):
            tbl = np.zeros(shape=refEntry.length,
                           dtype=self.LOCUS_SUMMARY_DTYPE)
            ra = tbl.view(np.recarray)
            ra.consensus[:] = list(noEvidenceConsensusCall(refEntry.sequence.tostring(),
                                                           options.noEvidenceConsensusCall))
            self.consensusByRefId[refId] = tbl

    def installInTable(self, locusSummary, tbl, rowNumber):
        # The consensus field in the table is type S8, so can only
        # store 8 characters.  On the odd chance that we detect an 9+
        # character insertion, we lose its first bases.
        tbl[rowNumber] = locusSummary._replace(consensus=locusSummary.consensus[-8:])[1:]


# We define both process and thread-based plurality collectors.
class PluralityResultCollectorProcess(PluralityResult, ResultCollectorProcess):
    pass

class PluralityResultCollectorThread(PluralityResult, ResultCollectorThread):
    pass


#
# Plugin API
#

# Pluggable module API for algorithms:
#  - Algorithm lives in a package
#  - Package must never fail to import, even if some of
#    its dependencies are not installed.
#  - Package must provide a main module exposing these top level
#    variables/methods:
#    - name                            = str
#    - availability                    = (bool, str)
#    - additionalDefaultOptions        = dict (string->value)
#    - compatibilityWithCmpH5(cmpH5)  -> (bool, str)
#    - slaveFactories                 -> bool -> (class, class)

__all__ = [ "name",
            "availability",
            "additionalDefaultOptions",
            "compatibilityWithCmpH5",
            "slaveFactories" ]

name = "Plurality"
availability = (True, "OK")

additionalDefaultOptions = { "referenceChunkOverlap"      : 0,
                             "variantCoverageThreshold"   : 3,
                             "variantConfidenceThreshold" : 20,
                             "coverage"                   : 250 }

def compatibilityWithCmpH5(cmpH5):
    return (True, "OK")

def slaveFactories(threaded):
    if threaded:
        return (PluralityWorkerThread,  PluralityResultCollectorThread)
    else:
        return (PluralityWorkerProcess, PluralityResultCollectorProcess)
