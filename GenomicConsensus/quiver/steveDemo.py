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
