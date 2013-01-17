import numpy as np, sys
from matplotlib.pyplot import *
from matplotlib.cm import binary

import ConsensusCore as cc
from GenomicConsensus.quiver import MutationScorer, VariantScorer
from GenomicConsensus.variants import Insertion, Deletion, Snv
from AlignmentHitStubs import AlignmentHitStub, _

referenceWindow = (1, 0, 100)

def makeScoreGrid(correctRead, incorrectRead, variant):
    """
    Build a grid of score(reference) - score(reference + variant)
    """
    correctScore = np.ndarray(shape=[20, 20])
    incorrectScore = np.ndarray(shape=[20, 20])

    for nCorrectReads in xrange(20):
        for nIncorrectReads in xrange(20):
            reads = nIncorrectReads*[incorrectRead] + nCorrectReads*[correctRead]
            vs = VariantScorer(referenceWindow)
            for read in reads:
                vs.addAlignment(read)
                correctScore[nCorrectReads, nIncorrectReads] = vs.scoreWithoutVariant(variant)
                incorrectScore[nCorrectReads, nIncorrectReads] = vs.scoreWithVariant(variant)
    return correctScore - incorrectScore


def demo1():
    print
    print "Demo part I: wherein we show that we have some sensitivity to overturn an "
    print "insertion call, based on the InsQV."

    insertionVariant = Insertion(1, 8,  8, "-", "A")

    #                                          012345678901234567
    AAread        = AlignmentHitStub(0, True, "GATCGATCAAGATCGATC",
                                              "GATCGATCAAGATCGATC",
                                      InsQV=_("                  "),
                                     SubsQV=_("                  "),
                                      DelQV=_("                  "),
                                     DelTag=_("NNNNNNNNNNNNNNNNNN"),
                                    MergeQV=_("                  "))

    #                                         01234567-8901234567
    AAAread     = AlignmentHitStub(0,  True, "GATCGATC-AAGATCGATC",
                                             "GATCGATCAAAGATCGATC",
                                     InsQV=_("        *          "),
                                    SubsQV=_("                   "),
                                     DelQV=_("                   "),
                                    DelTag=_("NNNNNNNNNNNNNNNNNNN"),
                                   MergeQV=_("                   "))

    Z = makeScoreGrid(correctRead=AAread,
                      incorrectRead=AAAread,
                      variant=insertionVariant)

    clf()
    imshow(Z>0, origin="lower", interpolation="nearest")
    ylabel("# reads with AA")
    xlabel("# reads with AAA")
    title("Insertion test: Quiver calls AA (red=TRUE)")


def demo2():
    print
    print "Demo part III: wherein we show that we have some sensitivity to overturn an "
    print "deletion call, based on the DelQV and DelTag."

    referenceWindow = (1, 0, 100)
    deletionVariant = Deletion(1, 8,  9, "A", "-")

    #                                         0123456789012345678
    AAread_correctDelTag = \
                   AlignmentHitStub(0, True, "GATCGATCAAAGATCGATC",
                                             "GATCGATC-AAGATCGATC",
                                     InsQV=_("                   "),
                                    SubsQV=_("                   "),
                                     DelQV=_("         *         "),
                                    DelTag=_("NNNNNNNNNANNNNNNNNN"),
                                   MergeQV=_("                   "))

    #                                         0123456789012345678
    AAread_incorrectDelTag = \
                   AlignmentHitStub(0, True, "GATCGATCAAAGATCGATC",
                                             "GATCGATC-AAGATCGATC",
                                     InsQV=_("                   "),
                                    SubsQV=_("                   "),
                                     DelQV=_("         *         "),
                                    DelTag=_("NNNNNNNNNNNNNNNNNNN"),
                                   MergeQV=_("                   "))


    #                                         0123456789012345678
    AAAread     = AlignmentHitStub(0,  True, "GATCGATCAAAGATCGATC",
                                             "GATCGATCAAAGATCGATC",
                                     InsQV=_("                   "),
                                    SubsQV=_("                   "),
                                     DelQV=_("                   "),
                                    DelTag=_("NNNNNNNNNNNNNNNNNNN"),
                                   MergeQV=_("                   "))

    Z_correctDelTag = makeScoreGrid(correctRead=AAAread,
                                    incorrectRead=AAread_correctDelTag,
                                    variant=deletionVariant)

    Z_incorrectDelTag = makeScoreGrid(correctRead=AAAread,
                                      incorrectRead=AAread_incorrectDelTag,
                                      variant=deletionVariant)
    clf()
    subplot(121)
    imshow(Z_correctDelTag>0, origin="lower", interpolation="nearest")
    ylabel("# reads with AAA")
    xlabel("# reads with AA")
    title("Correct DelTag")
    subplot(122)
    imshow(Z_incorrectDelTag>0, origin="lower", interpolation="nearest")
    ylabel("# reads with AAA")
    xlabel("# reads with AA")
    title("Incorrect DelTag")

    suptitle("Deletion tests: red=quiver calls AAA (correct)")


def demo3():
    print
    print "Demo part IV: wherein we show that we have some sensitivity to overturn an "
    print "SNV call, based on the SubsQV"

    referenceWindow = (1, 0, 100)
    substitutionVariant = Snv(1, 9,  10, "A", "T")

    #                                         0123456789012345678
    AAAread =      AlignmentHitStub(0, True, "GATCGATCAAAGATCGATC",
                                             "GATCGATCAAAGATCGATC",
                                     InsQV=_("                   "),
                                    SubsQV=_("                   "),
                                     DelQV=_("                   "),
                                    DelTag=_("NNNNNNNNNNNNNNNNNNN"),
                                   MergeQV=_("                   "))

    #                                         0123456789012345678
    ATAread =       AlignmentHitStub(0, True,"GATCGATCAAAGATCGATC",
                                             "GATCGATCATAGATCGATC",
                                     InsQV=_("                   "),
                                    SubsQV=_("         *         "),
                                     DelQV=_("                   "),
                                    DelTag=_("NNNNNNNNNNNNNNNNNNN"),
                                   MergeQV=_("                   "))

    Z = makeScoreGrid(correctRead=AAAread,
                      incorrectRead=ATAread,
                      variant=substitutionVariant)

    clf()
    subplot(111)
    imshow(Z>0, origin="lower", interpolation="nearest")
    ylabel("# reads with AAA")
    xlabel("# reads with ATA")
    title("Quiver calls AAA (red=TRUE)")


demo1()
raw_input()
demo2()
raw_input()
demo3()
