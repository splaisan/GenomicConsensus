from pbcore.io import CmpH5Reader, FastaReader
from GenomicConsensus.quiver.utils import *
from GenomicConsensus.quiver.model import *
import argparse, random, re, string
import ConsensusCore as cc


styleTbl = string.maketrans("MRID", "|*  ")

def styled(t):
    """
    Gusfield transcript -> exonerate+
    """
    return t.translate(styleTbl)

def windowFromString(s):
    if s == None:
        return None
    m = re.match("(.*):(.*)-(.*)", s)
    assert m
    refId    = int(m.group(1))
    refStart = int(m.group(2))
    refEnd   = int(m.group(3))
    return (refId, refStart, refEnd)

def getopts():
    parser = argparse.ArgumentParser(description="Demo program")
    parser.add_argument("--coverageDepth", "-D",
                        action="store",
                        type=int,
                        default=10)
    parser.add_argument("--referenceWindow", "-w",
                        action="store",
                        type=windowFromString,
                        default=(1,2000,2040))
    parser.add_argument("--referenceFasta", "-r",
                        action="store",
                        required=True,
                        type=str)
    parser.add_argument("inputFilename",
                        type=str,
                        help="The input cmp.h5 file")

    options = parser.parse_args()
    return options

def printGuess(truth, guess):
    aln = cc.Align(truth, guess)
    print "Truth:        " + aln.Target()
    print "(Comparison): " + styled(aln.Transcript())
    print "Guess:        " + aln.Query()


def printBanner(msg):
    print
    print msg
    print "-"*60

def main():
    options = getopts()
    c = CmpH5Reader(options.inputFilename)
    refString = next(iter(FastaReader(options.referenceFasta))).sequence
    PARAMS = AllQVsModel.trainedParams1().qvModelParams

    winId, winStart, winEnd = options.referenceWindow
    spanningFwdReads = [a for a in c.readsInRange(*options.referenceWindow)
                        if a.spansReferenceRange(winStart, winEnd)
                        if a.isForwardStrand]
    chosenReads = spanningFwdReads[:options.coverageDepth]
    clippedReads = [a.clippedTo(winStart, winEnd) for a in chosenReads]

    print
    print "True sequence:"
    trueTemplate = refString[winStart:winEnd]
    print trueTemplate
    print

    readStrings = [a.read(aligned=False, orientation="genomic")
                   for a in clippedReads]
    print
    print "Reads:"
    for s in readStrings:
        print s

    # POA is too good--only 1 error.  Let's sprinkle a couple errors
    # on the true template so we start from a suitably bad guess
    random.seed(4006)
    randomMutations = random.sample(list(uniqueSingleBaseMutations(trueTemplate)), 10)

    mutatedTemplate =  cc.ApplyMutations(randomMutations, trueTemplate)

    r = cc.SparseSseQvRecursor()
    mms = cc.SparseSseQvMultiReadMutationScorer(r, PARAMS, mutatedTemplate)
    for cr in clippedReads:
        mms.AddRead(AllQVsModel.extractFeatures(cr), int(cr.RCRefStrand))

    SEPARATION = 7
    NEIGHBORHOOD = 15
    favorableMutationsAndScores = []

    for round_ in range(10):
        printBanner("Round %d" % round_)
        printGuess(trueTemplate, mms.Template())

        mutationsToTry = uniqueSingleBaseMutations(mms.Template())
        favorableMutationsAndScores = \
            [(m, mms.Score(m)) for m in
             filter(mms.FastIsFavorable, mutationsToTry)]

        if not favorableMutationsAndScores:
            print
            print "Done!"
            print
            break

        else:
            bestMutationsAndScores = bestSubset(favorableMutationsAndScores, SEPARATION)
            bestMutations = map(fst, bestMutationsAndScores)

            print
            print "Applying mutations:"
            for m, s in bestMutationsAndScores:
                print m.ToString(), s

            mms.ApplyMutations(bestMutations)




if __name__ == '__main__':
    main()
