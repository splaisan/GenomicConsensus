from __future__ import absolute_import


import collections, h5py, math, logging, numpy as np, os, pprint
import sys, time
from .. import reference
from ..options import options
from ..Worker import WorkerProcess, WorkerThread
from ..ResultCollector import ResultCollectorProcess, ResultCollectorThread
from pbcore.io import rangeQueries, GffWriter
from ..io.fastx import FastaWriter, FastqWriter
from ..io.VariantsGffWriter import VariantsGffWriter

try:
    import ConsensusCore as cc
    from GenomicConsensus.quiver.utils import *
    from GenomicConsensus.quiver.model import *
    if cc.Version.VersionString() == "0.4.0":
        availability = (True, "OK")
    else:
        availability = (False, "Need ConsensusCore version 0.4.0")
except ImportError:
    availability = (False, "ConsensusCore not installed---required for Quiver algorithm")


#
# Tuning parameters, here for now
#
POA_COVERAGE    = 11
MUTATION_SEPARATION = 7

#
# The type used in the wire protocol.  Note that whereas plurality
# passes a packet per locus, we pass a bigger packet per reference
# window (1000bp, typically).
#
QuiverWindowSummary = collections.namedtuple("QuiverWindowSummary",
                                             ("referenceId",
                                              "referenceStart",
                                              "referenceEnd",
                                              "consensus",      # string
                                              "qv",             # numpy array of uint8
                                              "variants"))      # list of Variant

def domain(referenceWindow):
    """
    Calculate the "domain" within the referenceWindow---the region
    where we are allowed to call variants, since there is overlap
    between windows.
    """
    refId, winStart, winEnd = referenceWindow
    _, refGroupStart, refGroupEnd = \
        options.referenceWindow or \
        (refId, 0, reference.byId[refId].length)
    domainStart = winStart if winStart==refGroupStart else min(winStart + 5, winEnd)
    domainEnd   = winEnd   if winEnd  ==refGroupEnd   else max(winEnd   - 5, winStart)
    return (domainStart, domainEnd)


def extractClippedAlignments(referenceWindow, coverage, alnHits):
    """
    Extract up to `coverage` clipped alignments from alnHits, where
    spanning alignments are used first and then long nonspanning
    alignments.  The return value is a tuple of the lists,
      (clippedSpanningAlns, clippedNonSpanningAlns).
    """
    refId, refStart, refEnd = referenceWindow
    spanningAlns = [a for a in alnHits
                     if a.spansReferenceRange(refStart, refEnd)]
    # HACK below is to work around bug 21232 until it can be properly fixed
    nonSpanningAlns = [a for a in alnHits
                       if not a.spansReferenceRange(refStart, refEnd)
                       if referenceSpanWithinWindow(referenceWindow, a) > 0]  # <- HACK

    # We always prefer spanning reads, and among partial passes,
    # we favor those with a longer read length within the window.
    nonSpanningAlns.sort(key=lambda aln: -referenceSpanWithinWindow(referenceWindow, aln))

    def isStumpy(clippedAln):
        # Predicate for "stumps", where for example the aligner may
        # have inserted a large gap, such that while the alignment
        # technically spans the window, it may not have any read
        # content therein:
        #
        #   Ref   ATGATCCAGTTACTCCGATAAA
        #   Read  ATG---------------TA-A
        #   Win.     [              )
        #
        minimumReadLength = 0.1*(refEnd-refStart)
        return (clippedAln.readLength < minimumReadLength)

    # Some ugly logic for taking the best reads available until we
    # reach the desired coverage. Couldn't get this to fit into a nice
    # takewhile/filter functional style.
    clippedSpanningAlns = []
    clippedNonSpanningAlns = []
    for aln in spanningAlns:
        if len(clippedSpanningAlns) == coverage: break
        ca = aln.clippedTo(refStart, refEnd)
        if not isStumpy(ca):
            clippedSpanningAlns.append(ca)
    for aln in nonSpanningAlns:
        if len(clippedSpanningAlns) + len(clippedNonSpanningAlns) == coverage: break
        ca = aln.clippedTo(refStart, refEnd)
        clippedNonSpanningAlns.append(ca)

    return (clippedSpanningAlns, clippedNonSpanningAlns)

def variantsFromConsensus(refWindow, refSequenceInWindow, cssSequenceInWindow,
                          cssQvInWindow=None, siteCoverage=None, aligner="affine"):
    """
    Compare the consensus and the reference in this window, returning
    a list of variants.
    """
    refId, refStart, refEnd = refWindow

    if aligner == "affine":
        align = cc.AlignWithAffineGapPenalty
    else:
        align = cc.Align

    ga = align(refSequenceInWindow, cssSequenceInWindow)
    variants =  variantsFromAlignment(ga, refWindow)

    cssPosition = cc.TargetToQueryPositions(ga)

    for v in variants:
        # HACK ALERT: we are not really handling the scoring of
        # deletions at the end of the window correctly here.  The
        # correct fix will be to test an insertion at the end of the
        # template, lengthening the consensusQv by one.
        cssStart = min(cssPosition[v.refStart-refStart], len(cssQvInWindow)-1)

        if siteCoverage  != None: v.coverage   = siteCoverage[v.refStart-refStart]
        if cssQvInWindow != None: v.confidence = cssQvInWindow[cssStart]

    return variants


def dumpEvidence(evidenceDumpBaseDirectory,
                 refWindow, refSequenceInWindow,
                 clippedSpanningAlnsInWindow, clippedNonSpanningAlnsInWindow,
                 poaConsensusInWindow, quiverConsensusInWindow, scoresMatrix):
    # Format of evidence dump:
    # evidence_dump/
    #   ref000001/
    #     0-1005/
    #       reference.fa
    #       reads.fa
    #       poa-consensus.fa
    #       quiver-consensus.fa
    #       quiver-scores.h5
    #     995-2005/
    #       ...
    join = os.path.join
    refId, refStart, refEnd = refWindow
    refName = reference.idToName(refId)
    windowDirectory = join(evidenceDumpBaseDirectory,
                           refName,
                           "%d-%d" % (refStart, refEnd))
    logging.info("Dumping evidence to %s" % (windowDirectory,))

    if os.path.exists(windowDirectory):
        raise Exception, "Evidence dump does not expect directory %s to exist." % windowDirectory
    os.makedirs(windowDirectory)
    refFasta             = FastaWriter(join(windowDirectory, "reference.fa"))
    readsFasta           = FastaWriter(join(windowDirectory, "reads.fa"))
    poaConsensusFasta    = FastaWriter(join(windowDirectory, "poa-consensus.fa"))
    quiverConsensusFasta = FastaWriter(join(windowDirectory, "quiver-consensus.fa"))

    windowName = refName + (":%d-%d" % (refStart, refEnd))
    refFasta.writeRecord(windowName, refSequenceInWindow)
    refFasta.close()

    poaConsensusFasta.writeRecord(windowName + "|poa", poaConsensusInWindow)
    poaConsensusFasta.close()

    quiverConsensusFasta.writeRecord(windowName + "|quiver", quiverConsensusInWindow)
    quiverConsensusFasta.close()

    # TODO: add row/col names to make the matrix interpretable
    quiverScoreFile = h5py.File(join(windowDirectory, "quiver-scores.h5"))
    ds = quiverScoreFile.create_dataset("Scores", data=scoresMatrix)
    quiverScoreFile.close()

    for aln in (clippedSpanningAlnsInWindow + clippedNonSpanningAlnsInWindow):
        readsFasta.writeRecord(aln.readName, aln.read(orientation="genomic", aligned=False))
    readsFasta.close()

class QuiverWorker(object):

    def onStart(self):
        if options.parameters == "best":
            self.params = ParameterSet.bestAvailable(self._inCmpH5)
        else:
            self.params = ParameterSet.fromString(options.parameters)
        self.model  = self.params.model

        # If the parameter set was chosen automatically by "best" and
        # the model is not AllQVs, we should warn the user that this
        # could result in a degradation in performance.
        if options.parameters=="best" and self.model != AllQVsModel:
            logging.warn(
                "This .cmp.h5 file lacks some of the QV data tracks that are required " +
                "for optimal performance of the Quiver algorithm.  For optimal results" +
                " use the ResequencingQVs workflow in SMRTPortal with bas.h5 files "    +
                "from an instrument using software version 1.3.1 or later.")


    def onChunk(self, referenceWindow, alnHits):
        refId, refStart, refEnd = referenceWindow
        domainStart, domainEnd = domain(referenceWindow)
        domainLen = domainEnd - domainStart
        refSequence = reference.byId[refId].sequence[refStart:refEnd].tostring()
        domainRefSequence = reference.byId[refId].sequence[domainStart:domainEnd].tostring()

        clippedSpanningAlns, clippedNonSpanningAlns = \
            extractClippedAlignments(referenceWindow, options.coverage, alnHits)
        clippedAlns = clippedSpanningAlns + clippedNonSpanningAlns
        coverageReport = "[span=%d, nonspan=%d]" % \
            (len(clippedSpanningAlns), len(clippedNonSpanningAlns))
        siteCoverage = rangeQueries.getCoverageInRange(self._inCmpH5, referenceWindow,
                                                       rowNumbers=[a.rowNumber
                                                                   for a in clippedAlns])

        # In the case where there are no spanning reads, we really
        # cannot come up with a good consensus if the sequence is
        # highly divergent from the reference.  By default, the
        # behavior is to "nocall" the window in question, but if the
        # user sets --noEvidenceConsensusCall="reference", the
        # reference will be used instead.

        if len(clippedSpanningAlns) < 3:
            if options.noEvidenceConsensusCall == "nocall":
                domainCss = "N"*domainLen
                domainQv  = np.zeros(domainLen, dtype=np.uint8)
                domainVariants = []
                return QuiverWindowSummary(refId, refStart, refEnd,
                                           domainCss, domainQv, domainVariants)
            elif len(clippedNonSpanningAlns) == 0:
                domainCss = domainRefSequence
                domainQv  = np.zeros(domainLen, dtype=np.uint8)
                domainVariants = []
                return QuiverWindowSummary(refId, refStart, refEnd,
                                           domainCss, domainQv, domainVariants)
            # ... otherwise, we continue, using the reference as the
            # starting point for quiver's estimation.

        # Load the bits that POA cares about.
        forwardStrandSequences = [a.read(orientation="genomic", aligned=False)
                                  for a in clippedSpanningAlns]
        # Load the bits that quiver cares about
        mappedReads = [self.model.extractMappedRead(ca, refStart)
                       for ca in clippedAlns]

        # Beyond this point, no more reference will be made to the
        # pbcore alignment objects.

        # If there are any spanning reads, we can construct a POA
        # consensus as a starting point.  If there are no spanning
        # reads, we use the reference.
        if len(clippedSpanningAlns) == 0:
            css = refSequence
            numPoaVariants = None
        else:
            # We calculate a quick consensus estimate using the
            # Partial Order Aligner (POA), which typically gets us to
            # > 99% accuracy.  We have to lift the reference
            # coordinates into coordinates within the POA consensus.
            p = cc.PoaConsensus.FindConsensus(forwardStrandSequences[:POA_COVERAGE])
            ga = cc.Align(refSequence, p.Sequence())
            numPoaVariants = ga.Errors()
            css = poaCss = p.Sequence()

            # Lift reference coordinates onto the POA consensus coordinates
            queryPositions = cc.TargetToQueryPositions(ga)
            mappedReads = [lifted(queryPositions, mr) for mr in mappedReads]

        # Load the reads, including QVs, into a MutationScorer, which is a
        # principled and fast way to test potential refinements to our
        # consensus sequence.
        mms = cc.SparseSseQvMultiReadMutationScorer(self.params.quiverConfig, css)
        for mr in mappedReads:
            mms.AddRead(mr)

        # Test mutations, improving the consensus
        quiverCss, quiverConverged = refineConsensus(mms)
        css = quiverCss
        cssQv = consensusConfidence(mms)

        ga = cc.Align(refSequence, css)
        targetPositions = cc.TargetToQueryPositions(ga)

        # Calculate variants
        variants = variantsFromConsensus(referenceWindow, refSequence, css,
                                         cssQv, siteCoverage, options.aligner)

        numQuiverVariants = len(variants)

        # Restrict the consensus and variants to the domain.
        cssDomainStart = targetPositions[domainStart-refStart]
        cssDomainEnd   = targetPositions[domainEnd-refStart]
        domainCss = css[cssDomainStart:cssDomainEnd]
        domainQv = cssQv[cssDomainStart:cssDomainEnd]
        domainVariants = [ v for v in variants
                           if domainStart <= v.refStart < domainEnd ]

        poaReport    = "POA unavailable" if numPoaVariants==None \
                       else ("%d POA variants" % numPoaVariants)
        quiverReport = "%d Quiver variants" % numQuiverVariants

        logging.info("%s: %s, %s %s" %
                     (referenceWindow, poaReport, quiverReport, coverageReport))
        if not quiverConverged:
            logging.info("%s: Quiver did not converge to MLE" % (referenceWindow,))

        shouldDumpEvidence = \
            ((options.dumpEvidence == "all") or
             (options.dumpEvidence == "variants") and (numQuiverVariants > 0))

        if shouldDumpEvidence:
            quiverScores = scoreMatrix(mms)
            dumpEvidence(options.evidenceDirectory,
                         referenceWindow, refSequence,
                         clippedSpanningAlns, clippedNonSpanningAlns,
                         poaCss, quiverCss, quiverScores)

        return QuiverWindowSummary(refId, refStart, refEnd,
                                   domainCss, domainQv, domainVariants)


ContigConsensusChunk = collections.namedtuple("ContigConsensusChunk",
                                              ("refStart", "refEnd", "consensus", "qv"))

class QuiverResultCollector(object):

    def onStart(self):
        self.allVariants = []
        # this is a map of refId -> [ContigConsensusChunk];
        self.consensusChunks = collections.defaultdict(list)

    def onResult(self, result):
        assert result.consensus != None and isinstance(result.consensus, str)
        self.allVariants += result.variants
        cssChunk = ContigConsensusChunk(result.referenceStart,
                                        result.referenceEnd,
                                        result.consensus,
                                        result.qv)
        self.consensusChunks[result.referenceId].append(cssChunk)

    def onFinish(self):
        # 1. GFF output.
        if options.gffOutputFilename:
            # Dictionary of refId -> [Variant]
            filteredVariantsByRefId = collections.defaultdict(list)
            for v in self.allVariants:
                if (v.confidence > options.variantConfidenceThreshold and
                    v.coverage   > options.variantCoverageThreshold):
                    filteredVariantsByRefId[v.refId].append(v)
            self.writeVariantsGff(options.gffOutputFilename, filteredVariantsByRefId)

        # 2. FASTA output.  Currently unoptimized--will choke hard on
        # very large references.
        if options.fastaOutputFilename:
            writer = FastaWriter(options.fastaOutputFilename)
            for refId, unsortedChunks in self.consensusChunks.iteritems():
                chunks = sorted(unsortedChunks)
                css = "".join(chunk.consensus for chunk in chunks)
                quiverHeader = reference.idToHeader(refId) + "|quiver"
                writer.writeRecord(quiverHeader, css)

        # 3. FASTQ output
        if options.fastqOutputFilename:
            writer = FastqWriter(options.fastqOutputFilename)
            for refId, unsortedChunks in self.consensusChunks.iteritems():
                chunks = sorted(unsortedChunks)
                css = "".join(chunk.consensus for chunk in chunks)
                qv = np.concatenate([chunk.qv for chunk in chunks])
                quiverHeader = reference.idToHeader(refId) + "|quiver"
                writer.writeRecord(quiverHeader, css, qv)


    def writeVariantsGff(self, filename, filteredVariantsByRefId):
        writer = VariantsGffWriter(options.gffOutputFilename,
                                   " ".join(sys.argv),
                                   reference.byId.values())
        for id in reference.byId:
            for v in filteredVariantsByRefId[id]:
                writer.writeRecord(v.toGffRecord())


#
# Slave process/thread classes
#
class QuiverWorkerProcess(QuiverWorker, WorkerProcess): pass
class QuiverWorkerThread(QuiverWorker, WorkerThread): pass
class QuiverResultCollectorProcess(QuiverResultCollector, ResultCollectorProcess): pass
class QuiverResultCollectorThread(QuiverResultCollector, ResultCollectorThread):  pass


#
# Plugin API
#
__all__ = [ "name",
            "availability",
            "additionalDefaultOptions",
            "compatibilityWithCmpH5",
            "slaveFactories" ]

name = "Quiver"
additionalDefaultOptions = { "referenceChunkOverlap"      : 5,
                             "variantCoverageThreshold"   : 11,
                             "variantConfidenceThreshold" : 40,
                             "coverage"                   : 100,
                             "parameters"                 : "best" }

def compatibilityWithCmpH5(cmpH5):
    if cmpH5.readType != "standard":
        return (False, "The Quiver algorithm requires a cmp.h5 file containing standard (non-CCS) reads.")
    elif (options.parameters != "best" and
        not ParameterSet.fromString(options.parameters).model.isCompatibleWithCmpH5(cmpH5)):
        return (False, "This Quiver parameter set requires QV features not available in this .cmp.h5 file.")
    else:

        return (True, "OK")

def slaveFactories(threaded):
    # By default we use slave processes. The tuple ordering is important.
    if threaded:
        return (QuiverWorkerThread,  QuiverResultCollectorThread)
    else:
        return (QuiverWorkerProcess, QuiverResultCollectorProcess)
