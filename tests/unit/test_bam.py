import pbcore.data as D
from pbcore.io import CmpH5Reader
from GenomicConsensus.io import PacBioBamReader, BamReader

from numpy.testing import (assert_array_equal        as ARRAY_EQ,
                           assert_array_almost_equal as ARRAY_SIM)
from nose.tools import assert_equal as EQ, assert_raises
import numpy as np


class TestBam(object):

    def __init__(self):
        bamFname, cmpFname = D.getBamAndCmpH5()
        lambdaFasta = D.getLambdaFasta()

        self.b = PacBioBamReader(bamFname, lambdaFasta)
        self.c = CmpH5Reader(cmpFname)
        self.bBasic = BamReader(bamFname)

        # Note that sorting orders are not generally the same... BAM
        # sorts + alns before - alns, when there is a tie on tStart;
        # we don't do this in cmp.h5 (we next sort on tEnd).  However
        # in this file there are no ties on tStart.
        self.bAlns = list(self.b)
        self.bFwd = self.bAlns[0]
        self.bRev = self.bAlns[1]

        self.cAlns = list(self.c)
        self.cFwd = self.cAlns[0]
        self.cRev = self.cAlns[1]

        self.cFwdClipped = self.cFwd.clippedTo(10, 60)
        self.bFwdClipped = self.bFwd.clippedTo(10, 60)
        self.cRevClipped = self.cRev.clippedTo(310, 360)
        self.bRevClipped = self.bRev.clippedTo(310, 360)



    @staticmethod
    def compareAlns(b, c):
        EQ(c.readName        , b.readName)
        EQ(c.isForwardStrand , b.isForwardStrand)
        EQ(c.isReverseStrand , b.isReverseStrand)
        EQ(c.rStart          , b.rStart)
        EQ(c.rEnd            , b.rEnd)
        EQ(c.tStart          , b.tStart)
        EQ(c.tEnd            , b.tEnd)
        EQ(c.HoleNumber      , b.HoleNumber)

        for o in ["native", "genomic"]:

            for s in ("cigar", "exonerate", "exonerate+", "gusfield"):
                EQ(c.transcript(style=s),
                   b.transcript(style=s))

            ARRAY_EQ(c.referencePositions(orientation=o),
                     b.referencePositions(orientation=o))

            ARRAY_EQ(c.readPositions(orientation=o),
                     b.readPositions(orientation=o))

            for a in [True, False]:
                ARRAY_EQ(c.read(a, o),
                         b.read(a, o))

                ARRAY_EQ(c.reference(a, o),
                         b.reference(a, o))

                # MergeQV is actually not identical presently, due to the
                # clipping at 93 used for ASCII encoding in BAM
                for featureName in [ "DeletionQV",
                                     "InsertionQV",
                                     "SubstitutionQV" ]:
                    ARRAY_EQ(c.pulseFeature(featureName, a, o),
                             b.pulseFeature(featureName, a, o))


        # DelTag in "genomic" orientation is arguably done incorrectly by the
        # cmp.h5 reader (it doesn't RC the bases), so we don't compare that.
        ARRAY_EQ(c.pulseFeature("DeletionTag", aligned=False, orientation="native"),
                 b.pulseFeature("DeletionTag", aligned=False, orientation="native"))



    def testForwardStrandAln(self):
        self.compareAlns(self.cFwd, self.bFwd)


    def testReverseStrandAln(self):
        self.compareAlns(self.cRev, self.bRev)

    def testClipping(self):
        self.compareAlns(self.cFwdClipped, self.bFwdClipped)
        self.compareAlns(self.cRevClipped, self.bRevClipped)

    def testIndex(self):
        ARRAY_EQ([aln.tStart for aln in self.bAlns], self.b.pbi.tStart)
        ARRAY_EQ([aln.tEnd   for aln in self.bAlns], self.b.pbi.tEnd)
        ARRAY_EQ([aln.rStart for aln in self.bAlns], self.b.pbi.rStart)
        ARRAY_EQ([aln.rEnd   for aln in self.bAlns], self.b.pbi.rEnd)
        ARRAY_EQ([aln.MapQV  for aln in self.bAlns], self.b.pbi.MapQV)

    # def testClippingRegression(self):
    #     # Test a corner case
    #     clipC1 = self.cFwd.clippedTo(56, 58)
    #     clipB1 = self.bFwd.clippedTo(56, 58)
    #     self.compareAlns(clipC1, clipB1)

    # def testClippingExhaustively(self):
    #     pass

    def testIncorrectReference(self):
        bamFname, _ = D.getBamAndCmpH5()
        incorrectFasta = D.getTinyFasta()
        with assert_raises(Exception):
            f = BamReader(bamFname, incorrectFasta)

    # def testRangeQueries1(self):
    #     win = (0, 1000, 2000)
    #     expectedReadNames =

    # def testRangeQueries2(self):



BT = TestBam()
#BT.testForwardStrandAln()
#BT.testReverseStrandAln()
