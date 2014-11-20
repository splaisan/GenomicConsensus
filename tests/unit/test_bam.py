import pbcore.data as D
from pbcore.io import CmpH5Reader
from GenomicConsensus.io import BamReader

from numpy.testing import (assert_array_equal        as ARRAY_EQ,
                           assert_array_almost_equal as ARRAY_SIM)
from nose.tools import assert_equal as EQ

import numpy as np


class TestBam(object):

    def __init__(self):
        bamFname, cmpFname = D.getBamAndCmpH5()
        self.b = BamReader(bamFname)
        self.c = CmpH5Reader(cmpFname)

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

        for orientation in ["native", "genomic"]:
            for aligned in [True, False]:
                ARRAY_EQ(c.read(aligned, orientation),
                         b.read(aligned, orientation))

        # MergeQV is actually not identical presently, due to the
        # clipping at 93 used for ASCII encoding in BAM
        for featureName in [ "DeletionQV",
                             "InsertionQV",
                             "SubstitutionQV",
                             "DeletionTag" ]:

            for orientation in ["native", "genomic"]:
                for aligned in [True, False]:
                    ARRAY_EQ(c.pulseFeature(featureName, aligned, orientation),
                             b.pulseFeature(featureName, aligned, orientation))

        ARRAY_EQ(c.pulseFeature("DeletionTag", aligned=True, orientation="genomic"),
                 b.pulseFeature("DeletionTag", aligned=True, orientation="genomic"))

        for orientation in ["genomic", "native"]:
            ARRAY_EQ(c.referencePositions(orientation),
                     b.referencePositions(orientation))


    def testForwardStrandAln(self):
        self.compareAlns(self.cFwd, self.bFwd)


    def testReverseStrandAln(self):
        self.compareAlns(self.cRev, self.bRev)


    def testClipping(self):
        pass


BT = TestBam()
#BT.testForwardStrandAln()
#BT.testReverseStrandAln()


x = BT.cRev.DeletionQV(orientation="native", aligned=True)
y = BT.bRev.DeletionQV(orientation="native", aligned=True)


x = BT.cRev.pulseFeature("DeletionTag", orientation="genomic", aligned=True)
y = BT.bRev.pulseFeature("DeletionTag", orientation="genomic", aligned=True)


ARRAY_EQ(BT.cRev.read(True, "native"),
         BT.bRev.read(True, "native"))
