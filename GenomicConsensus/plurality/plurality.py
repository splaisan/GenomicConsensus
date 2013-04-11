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
from itertools import izip
from collections import Counter
from ..utils import *
from .. import reference
from ..options import options
from ..Worker import WorkerProcess, WorkerThread
from ..ResultCollector import ResultCollectorProcess, ResultCollectorThread
from ..consensus import *
from ..variants import *


#
# --------------- Configuration ----------------------
#

class PluralityConfig(object):
    """
    Plurality configuration options
    """
    def __init__(self,
                 minMapQV=10,
                 minCoverage=3,
                 minConfidence=40,
                 noEvidenceConsensus="nocall"):
        self.minMapQV            = minMapQV
        self.minCoverage         = minCoverage
        self.minConfidence       = minConfidence
        self.noEvidenceConsensus = noEvidenceConsensus
        self.realignHomopolymers = False # not available yet


#
# -----------  The actual algorithm code -------------
#

def pluralityConsensusAndVariants(refWindow, referenceSequenceInWindow, alns,
                                  pluralityConfig):

    """
    Compute (Consensus, [Variant]) for this window, using the given
    `alns`, by applying a straightforward column-oriented consensus
    calling algorithm.

    If the consensus cannot be called for a base, "N" will be placed
    in the consensus sequence for that position.

    If `realignHomopolymers` is True, alignment gaps will be shuffled
    in homopolymer regions in an attempt to maximize variant detection
    sensitivity (not yet implemented, and may never be).
    """
    _, refStart, refEnd = refWindow
    windowSize = refEnd - refStart
    assert len(referenceSequenceInWindow) == windowSize

    #
    # Build up these arrays in reference coordinates.
    #
    consensusSequence_   = []
    consensusFrequency_  = []
    effectiveCoverage_   = []
    noEvidenceConsensusFactory = \
        noEvidenceConsensusFactoryByName[pluralityConfig.noEvidenceConsensus]
    noEvidenceConsensus = noEvidenceConsensusFactory(refWindow, referenceSequenceInWindow)

    baseCallsMatrix = tabulateBaseCalls(refWindow, alns)

    for j in xrange(0, windowSize):
        counter = Counter(baseCallsMatrix[:, j])
        if "" in counter: counter.pop("")

        siteEffectiveCoverage = sum(counter.itervalues())
        if siteEffectiveCoverage < pluralityConfig.minCoverage:
            siteConsensusFrequency = siteEffectiveCoverage
            siteConsensusSequence  = noEvidenceConsensus.sequence[j]
        else:
            siteConsensusFrequency, siteConsensusSequence = max(zip(counter.values(),
                                                                    counter.keys()))
        # Replace explicit gaps with empty string
        if siteConsensusSequence == "-":
            siteConsensusSequence = ""

        consensusSequence_.append(siteConsensusSequence)
        consensusFrequency_.append(siteConsensusFrequency)
        effectiveCoverage_.append(siteEffectiveCoverage)

    consensusConfidence_ = computePluralityConfidence(consensusFrequency_,
                                                      effectiveCoverage_)

    #
    # Now we need to put everything in consensus coordinates
    #
    consensusLens = map(len, consensusSequence_)
    consensusSequence = "".join(consensusSequence_)
    consensusConfidence = np.repeat(consensusConfidence_, consensusLens)
    css = Consensus(refWindow, consensusSequence, consensusConfidence)

    #
    # Derive variants from reference-coordinates consensus
    #
    variants_ = _computeVariants(refWindow,
                                 referenceSequenceInWindow,
                                 consensusSequence_,
                                 consensusConfidence_,
                                 effectiveCoverage_,
                                 consensusFrequency_)
    variants = filterVariants(pluralityConfig.minCoverage,
                              pluralityConfig.minConfidence,
                              variants_)
    return (css, variants)


def _computeVariants(refWindow,
                     refSequenceInWindow,
                     consensusArray,
                     consensusConfidenceArray,
                     consensusCoverageArray,
                     consensusFrequencyArray,
                     minCoverage=3,
                     minConfidence=20):
    refId, refStart, refEnd = refWindow
    windowSize = refEnd - refStart
    assert len(refSequenceInWindow) == windowSize
    assert len(consensusArray) == windowSize

    vars = []
    for j in xrange(windowSize):
        refPos = j + refStart
        refBase = refSequenceInWindow[j]
        cssBases = consensusArray[j]
        conf = consensusConfidenceArray[j]
        cov  = consensusCoverageArray[j]
        freq = consensusFrequencyArray[j]

        if (refBase != cssBases):
            if cssBases == "":
                vars.append(Deletion(refId, refPos, refPos+1,
                                     refBase, cssBases, cov, conf, freq))
            else:
                if len(cssBases) > 1:
                    vars.append(Insertion(refId, refPos, refPos,
                                          "", cssBases[:-1], cov, conf, freq))
                if cssBases[-1] != refBase:
                    vars.append(Substitution(refId, refPos, refPos+1,
                                             refBase, cssBases[-1], cov, conf, freq))
    return sorted(vars)

def tabulateBaseCalls(refWindow, alns, realignHomopolymers=False):
    """
    Go through the reads and build up the structured baseCallsMatrix
    table, which tabulates the read bases occurring at each reference
    coordinate in each read.  This code is somewhat tricky, read carefully.
    """
    _, refStart, refEnd = refWindow
    windowSize = refEnd - refStart

    baseCallsMatrix = np.zeros(shape=(len(alns), windowSize), dtype="S8")

    for i, aln in enumerate(alns):
        aln = aln.clippedTo(refStart, refEnd)
        alnRef    = aln.reference(orientation="genomic")
        alnRead   = aln.read(orientation="genomic")
        if realignHomopolymers:
            alnRef, alnRead =  normalizeHomopolymerGaps(alnRef, alnRead)

        # Idea: scan through the ref, read; for each non-gap character
        # in ref, record all non-gap characters seen in read since
        # last ref gap.
        readBases = []
        accum = []
        for (refBase, readBase) in izip(alnRef, alnRead):
            if readBase != "-":
                readBases.append(readBase)
            if refBase != "-":
                basesForRefPos = "".join(readBases) if readBases else "-"
                accum.append(basesForRefPos)
                readBases = []
        s, e = (aln.referenceStart - refStart,
                aln.referenceEnd   - refStart)
        baseCallsMatrix[i, s:e] = accum
    return baseCallsMatrix


def computePluralityConfidence(consensusFrequency, effectiveCoverage):
    """
    Come up with a new simpler scheme.  Should only need frequency and
    coverage.
    """
    confidence = np.empty_like(consensusFrequency)
    confidence.fill(40)
    return confidence


#
# --------------  Plurality Worker class --------------------
#

class PluralityWorker(object):

    @property
    def pluralityConfig(self):
        return self._algorithmConfig

    def onStart(self):
        random.seed(42)

    def onChunk(self, workChunk):
        referenceWindow = workChunk.window
        noCallFn = noEvidenceConsensusFactoryByName[options.noEvidenceConsensusCall]
        refSeqInWindow = reference.sequenceInWindow(referenceWindow)

        if not workChunk.hasCoverage:
            noCallCss = noCallFn(referenceWindow, refSeqInWindow)
            return (referenceWindow, (noCallCss, []))

        rowNumbers = readsInWindow(self._inCmpH5, referenceWindow,
                                   depthLimit=options.coverage,
                                   minMapQV=options.minMapQV,
                                   strategy="longest",
                                   stratum=options.readStratum)
        alnHits = self._inCmpH5[rowNumbers]
        return (referenceWindow,
                pluralityConsensusAndVariants(referenceWindow, refSeqInWindow,
                                              alnHits, self.pluralityConfig))

# define both process and thread-based plurality callers
class PluralityWorkerProcess(PluralityWorker, WorkerProcess): pass
class PluralityWorkerThread(PluralityWorker, WorkerThread): pass

#
#  --------------------- Plugin API --------------------------------
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
#    - configure                      -> options -> cmph5 -> algorithm specific config object;
#                                        (can raise IncompatibleDataException)
#    - slaveFactories                 -> bool -> (class, class)

__all__ = [ "name",
            "availability",
            "additionalDefaultOptions",
            "configure",
            "slaveFactories" ]

name = "Plurality"
availability = (True, "OK")

additionalDefaultOptions = { "referenceChunkOverlap"      : 0 }


def slaveFactories(threaded):
    if threaded:
        return (PluralityWorkerThread,  ResultCollectorThread)
    else:
        return (PluralityWorkerProcess, ResultCollectorProcess)

def configure(options, cmpH5):
    pluralityConfig = PluralityConfig(minMapQV=options.minMapQV,
                                      minCoverage=options.minCoverage,
                                      minConfidence=options.minConfidence,
                                      noEvidenceConsensus=options.noEvidenceConsensusCall)
    return pluralityConfig
