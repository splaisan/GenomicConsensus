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

import numpy as np, itertools, logging, re
from collections import Counter

from GenomicConsensus.variants import *
from GenomicConsensus.utils import error_probability_to_qv
import ConsensusCore as cc

# Some lisp functions we want
fst   = lambda t: t[0]
snd   = lambda t: t[1]
third = lambda t: t[2]

def uniqueSingleBaseMutations(templateSequence, positions=None):
    """
    Return an iterator over all single-base mutations of a
    templateSequence that result in unique mutated sequences.
    """
    allBases = "ACGT"
    prevTplBase = None
    positions = positions or xrange(0, len(templateSequence))
    for tplStart in positions:
        tplBase = templateSequence[tplStart]
        # snvs
        for subsBase in allBases:
            if subsBase != tplBase:
                yield cc.Mutation(cc.SUBSTITUTION, tplStart, subsBase)
        # Insertions---only allowing insertions that are not cognate
        # with the previous base.
        for insBase in allBases:
            if insBase != prevTplBase:
                yield cc.Mutation(cc.INSERTION, tplStart, insBase)
        # Deletion--only allowed if refBase does not match previous tpl base
        if tplBase != prevTplBase:
            yield cc.Mutation(cc.DELETION, tplStart, "-")
        prevTplBase = tplBase

def nearbyMutations(mutations, tpl, neighborhoodSize):
    """
    Return mutations nearby the previously-tried mutations
    """
    mutationPositions = map(cc.Mutation.Position, mutations)
    nearbyPositions = set()
    for mp in mutationPositions:
        nearbyPositions.update(range(max(0, mp - neighborhoodSize),
                                     min(len(tpl), mp + neighborhoodSize)))
    return uniqueSingleBaseMutations(tpl, sorted(nearbyPositions))

def asFloatFeature(arr):
    return cc.FloatFeature(np.array(arr, dtype=np.float32))

def bestSubset(mutationsAndScores, separation):
    """
    Given a list of (mutation, score) tuples, this utility method
    greedily chooses the highest scoring well-separated elements.  We
    use this to avoid applying adjacent high scoring mutations, which
    are the rule, not the exception.  We only apply the best scoring one
    in each neighborhood, and then revisit the neighborhoods after
    applying the mutations.
    """
    input = mutationsAndScores[:]
    output = []

    while input:
        best = max(input, key=snd)
        output.append(best)
        nStart = best[0].Position() - separation
        nEnd   = best[0].Position() + separation
        for t in input[:]:
            if nStart <= t[0].Position() <= nEnd:
                input.remove(t)

    return output

def refineConsensus(mms, maxRounds=20):
    """
    Given a MultiReadMutationScorer, identify and apply favorable
    template mutations.  Return (consensus, didConverge) :: (str, bool)
    """
    SEPARATION = 10
    NEIGHBORHOOD = 15
    favorableMutationsAndScores = None

    converged = False
    for round_ in range(1, maxRounds):

        if favorableMutationsAndScores == None:
            mutationsToTry = uniqueSingleBaseMutations(mms.Template())
        else:
            favorableMutations = map(fst, favorableMutationsAndScores)
            mutationsToTry = nearbyMutations(favorableMutations, mms.Template(), NEIGHBORHOOD)

        favorableMutationsAndScores = \
            [(m, mms.Score(m)) for m in
             filter(mms.FastIsFavorable, mutationsToTry)]

        if favorableMutationsAndScores:
            bestMutations = map(fst, bestSubset(favorableMutationsAndScores, SEPARATION))
            logging.debug("Applying mutations: %s" % [mut.ToString() for mut in bestMutations])
            mms.ApplyMutations(bestMutations)
        else:
            # If we can't find any favorable mutations, our work is done.
            converged = True
            break

    logging.debug("Quiver: %d rounds" % round_)
    return (mms.Template(), converged)

def findDinucleotideRepeats(s):
    """
    string -> list( (length-2 string, (start_position, end_position)) )
    """
    repeatsFound = []
    allDinucs = [ a + b for a in "ACGT" for b in "ACGT" if a != b ]
    for dinuc in allDinucs:
        matchIter = re.finditer("(?:%s){3,}" % dinuc, s)
        for m in matchIter:
            repeatsFound.append((dinuc, m.span()))
    return repeatsFound

def refineDinucleotideRepeats(mms):
    """
    We have observed a couple instances where we call the consensus to
    be off the truth by +/- 1 dinucleotide repeat---we are getting
    trapped in an inferor local optimum, like so:

                           likelihood
    truth       ATATATAT      100
    quiver      AT--ATAT       90
    quiver+A    ATA-ATAT       85
    quiver+T    AT-TATAT       85

    To resolve this issue, we need to explore the likelihood change
    for wobbling on every dinucleotide repeat in the window.

    This code is not to be emulated!  It is a stopgap until we can
    bake efficient multi-base mutation testing into ConsensusCore.
    """
    def wobble(mms, dinuc, startPos, repeatCountDiff):
        b1, b2 = dinuc
        if repeatCountDiff == -1:
            mms.ApplyMutations([cc.Mutation(cc.DELETION, start, "-")])
            mms.ApplyMutations([cc.Mutation(cc.DELETION, start, "-")])
        elif repeatCountDiff == 1:
            mms.ApplyMutations([cc.Mutation(cc.INSERTION, start, b2)])
            mms.ApplyMutations([cc.Mutation(cc.INSERTION, start, b1)])

    for (dinuc, (start, end)) in findDinucleotideRepeats(mms.Template()):
        currentScore = mms.BaselineScore()
        wobble(mms, dinuc, start, 1)
        scoreP1 = mms.BaselineScore()

        wobble(mms, dinuc, start, -1)
        wobble(mms, dinuc, start, -1)
        scoreN1 = mms.BaselineScore()

        bestIdx = np.argmax([scoreN1, currentScore, scoreP1])

        if bestIdx == 0:
            continue
        elif bestIdx == 1:
            wobble(mms, dinuc, start, 1)
            continue
        else:
            wobble(mms, dinuc, start, 1)
            wobble(mms, dinuc, start, 1)
            continue

def consensusConfidence(mms, positions=None):
    """
    Returns an array of QV values reflecting the consensus confidence
    at each position specified.  If the `positions` argument is
    omitted, confidence values are returned for all positions in the
    consensus (mms.Template()).
    """
    # TODO: We should be using all mutations here, not just all
    # mutations leading to unique templates.  This should be a trivial
    # fix, post-1.4.
    css = mms.Template()
    allMutations = uniqueSingleBaseMutations(css, positions)
    cssQv = []

    for pos, muts in itertools.groupby(allMutations, lambda m: m.Position()):
        # Current score is '0'; exp(0) = 1
        altScores = [mms.FastScore(m) for m in muts]
        with np.errstate(over="ignore", under="ignore"):
            errProb = 1. - 1. / (1. + sum(np.exp(altScores)))
        cssQv.append(error_probability_to_qv(errProb))

    return np.array(cssQv, dtype=np.uint8)


def inverseMutations(windowStart, variant):
    """
    Given a (potentially multibase) variant, return the list of single
    base mutations that 'undo' the variant.
    """
    start = variant.refStart - windowStart
    length = len(variant)
    if isinstance(variant, Insertion):
        ms = [cc.Mutation(cc.DELETION, pos, "-")
              for pos in xrange(start, start+length)]
    elif isinstance(variant, Deletion):
        ms = [cc.Mutation(cc.INSERTION, start, base)
              for base in variant.refSequence]
    elif isinstance(variant, Substitution):
        ms = [cc.Mutation(cc.SUBSTITUTION, pos, base)
              for (pos, base) in zip(range(start, start+length),
                                     variant.refSequence)]
    else:
        raise Exception, "Should not reach here"
    return ms

def variantsFromAlignment(a, refWindow):
    """
    Extract the variants implied by a pairwise alignment to the
    reference.
    """
    variants = []
    refId, refStart, _ = refWindow
    refPos = refStart
    tbl = zip(a.Transcript(),
              a.Target(),
              a.Query())

    # We don't call variants where either the reference or css is 'N'
    grouper = lambda row: "N" if (row[1]=="N" or row[2]=="N") else row[0]
    runs = itertools.groupby(tbl, grouper)

    for code, run in runs:
        assert code in "RIDMN"
        run = list(run)
        ref = "".join(map(snd, run))
        refLen = len(ref) - Counter(ref)["-"]
        read = "".join(map(third, run))

        if code == "M" or code == "N":
            pass
        elif code == "R":
            assert len(read)==len(ref)
            variants.append(Substitution(refId, refPos, refPos+len(read), ref, read))
        elif code == "I":
            variants.append(Insertion(refId, refPos, refPos, "", read))
        elif code == "D":
            variants.append(Deletion(refId, refPos, refPos + len(ref), ref, ""))

        refPos += refLen

    return variants

def referenceSpanWithinWindow(referenceWindow, aln):
    """
    Helper function for sorting reads by their reference span
    after restriction to a window.
    """
    _, winStart, winEnd = referenceWindow
    return min(winEnd, aln.referenceEnd) - \
           max(winStart, aln.referenceStart)

def lifted(queryPositions, mappedRead):
    """
    Lift a mappedRead into a new coordinate system by using the
    position translation table `queryPositions`
    """
    newStart = queryPositions[mappedRead.TemplateStart]
    newEnd   = queryPositions[mappedRead.TemplateEnd]
    return cc.MappedRead(mappedRead.Features,
                         mappedRead.Strand,
                         newStart,
                         newEnd)

def scoreMatrix(mms):
    """
    Produce a matrix S where S_{ij} represents the score delta of
    mutation j against read i.

    TODO: add row/column names
    """
    css = mms.Template()
    allMutations = sorted(uniqueSingleBaseMutations(css))
    shape = (mms.NumReads(), len(allMutations))
    scoreMatrix = np.zeros(shape)

    for j, mut in enumerate(allMutations):
        mutScores = mms.Scores(mut)
        scoreMatrix[:, j] = mutScores

    return scoreMatrix

def domain(referenceWindow):
    """
    Calculate the "domain" within the referenceWindow---the region
    where we are allowed to call variants, since there is overlap
    between windows.
    """
    from GenomicConsensus.options import options
    from GenomicConsensus import reference

    refId, winStart, winEnd = referenceWindow
    _, refGroupStart, refGroupEnd = \
        options.referenceWindow or \
        (refId, 0, reference.byId[refId].length)
    domainStart = winStart if winStart==refGroupStart else min(winStart + 5, winEnd)
    domainEnd   = winEnd   if winEnd  ==refGroupEnd   else max(winEnd   - 5, winStart)
    return (domainStart, domainEnd)


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
