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

# Author: Jim Drake

# Rare variant caller algorithm.  Based on plurality caller.

from __future__ import absolute_import

import math, logging, numpy as np
from collections import defaultdict, namedtuple, OrderedDict
from bisect import bisect_left, bisect_right
from itertools import izip
from ..utils import error_probability_to_qv
from ..options import options
from .. import (io,
                reference)
from ..Worker import WorkerProcess, WorkerThread
from ..ResultCollector import ResultCollectorProcess, ResultCollectorThread
from ..alignment import AlignmentColumn
from ..plurality.plurality import PluralityResult, PluralityLocusSummary

try:
    from scipy.stats import binom
    availability = (True, "OK")
except:
    availability = (False, "Cannot import SciPy (required for rare variant analysis)")


# NOTES:
# This needs some rethinking.  We're returning multiple variants for a given
# column.   Each variant will have a frequency and confidence value associated
# with it.  The locus summary should reflect that somehow.  There may also
# need to be some re-design work here to better reflect the data structure
# needs of the rare variant algorithm.

# We're only using CCS reads here for a first pass.  To support a new downstream
# feature, correlated mutations, we may want to output read ids associated
# with each particular variant.
class RareAlignmentColumn(AlignmentColumn):
    def __init__(self, *args):
        super(RareAlignmentColumn,self).__init__(*args)

    def confidence(self):
        pass

    # NOTE: This differs from plurality.consensus in that it returns a list of locus
    # summaries instead of just one.  I'm leaning in favor of renaming this method
    # `summary` instead of `consensus`, which seems to make more sense.  However,
    # all the downstream processes look for a `consensus` field, so we'll leave
    # it for now.
    # This is really just a looser filter than plurality.
    # Isolate the test from the data structures, that way we can swap algos.
    def consensus(self):
        """
        Generates coverage information only for positions that have a variant.
        Binomial test runs using scipy.stats.binom.sf
        x = # of variant observations
        n = Total number of observations (coverage)
        f = x / n
        p = Probability of success, or that the variant observed is a sequencing error.
            For CCS reads this is presumed to be 0.01
        for each snippit:
            if n > 500 and f > 0.01:
                pval = binom.sf(x, n, p=0.01)
                save to output, calculate confidence based on pval
        """
        coverage = self.coverage()
        loci = []

        # Don't bother processing if they all match the reference.
        saf = self._snippetsAndFrequencies
        if self.referenceBase in saf and saf[self.referenceBase] == coverage:
            return loci

        # Rare variants are 1% < freq < 50% at a coverage > 500. Dominant
        # alleles (i.e., 50% < freq) are also detected and reported.
        for snippet, count in self.orderedSnippetsAndFrequencies:
            freq = count / float(coverage)
            if coverage > options.variantCoverageThreshold and freq > 0.01:
                # TODO: modify to use QV values (defaults to nominal CCS error rate)
                pval = binom.sf(count, coverage, 0.01)
                loci.append(PluralityLocusSummary(
                                self.referenceId,
                                self.referencePos,
                                coverage,
                                snippet,
                                error_probability_to_qv(pval),
                                count))

        # TODO: needed?  only using in unit tests right now
        # tie-breaking sort, same as alignment.AlignmentColumn
        loci.sort(key=lambda pls: (pls.consensusFrequency, pls.consensus))
        return loci

class RareCaller(object):

    # The type of the result that will be passed on to the result
    # collector is:
    #     (ReferenceWindow, list of (Locus, Tuple)]
    # where each Tuple must be convertible to a numpy array with dtype
    # as below:

    def onChunk(self, referenceWindow, alnHits):
        return (referenceWindow, self.rare(referenceWindow, alnHits))

    @staticmethod
    def rare(referenceWindow, alignments):
        """
        Modeled after the PluralityCaller.plurality method

        Input: referenceWindow, iterable of alignmentHits
        Output: list of (Locus, PluralityLocusSummary)

        Like plurality, we compute a tally of snippits and frequencies.
        Then perform a binomial analysis of the table and return a list of
        variants based on cutoffs on both the frequency and p-value.
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
                        alignmentColumns[refId, refPos] = RareAlignmentColumn(refId, refPos, refBase)
                    alignmentColumns[refId, refPos].addReadSnippet("".join(readBases))
                    readBases = []


        return [(locus, consensus)
                for locus, algCol in alignmentColumns.iteritems()
                for consensus in algCol.consensus() ]

class RareResult(PluralityResult):
    """
    We override some data structures and methods here in order to support the
    data coming from rare variants.  Consumers downstream need to be aware of
    this.  Rare variants will only create a GFF file, so only the GFF consumer
    needs to be aware of the changes.
    """
    def initTable(self, refId, refEntry):
        self.consensusByRefId[refId] = dict()

    def installInTable(self, locusSummary, tbl, rowNumber):
        # TODO: Test memory footprint
        if rowNumber not in tbl:
            tbl[rowNumber] = list()
        # filter out indels here
        if len(locusSummary.consensus) == 1:
            tbl[rowNumber].append(locusSummary)


# define both process and thread-based rare variant callers
class RareWorkerProcess(RareCaller,WorkerProcess):
    pass

class RareWorkerThread(RareCaller,WorkerThread):
    pass

# define both process and thread-based plurality collectors
class RareResultCollectorProcess(RareResult, ResultCollectorProcess):
    pass

class RareResultCollectorThread(RareResult, ResultCollectorThread):
    pass


#
# Plugin API
#
__all__ = [ "name",
            "availability",
            "additionalDefaultOptions",
            "compatibilityWithCmpH5",
            "slaveFactories" ]

name = "Rare variant analysis"

additionalDefaultOptions = { "referenceChunkOverlap"      : 0,
                             "variantCoverageThreshold"   : 500,
                             "variantConfidenceThreshold" : 20 }

def compatibilityWithCmpH5(cmpH5):
    # TODO: check whether the cmp.h5 is CCS
    return (True, "OK")

def slaveFactories(threaded):
    if threaded:
        return (RareWorkerThread,  RareResultCollectorThread)
    else:
        return (RareWorkerProcess, RareResultCollectorProcess)
