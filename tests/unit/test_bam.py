import pbcore.data as D
from pbcore.io import CmpH5Reader
from GenomicConsensus.io import BamReader

from numpy.testing import (assert_array_equal        as ARRAY_EQ,
                           assert_array_almost_equal as ARRAY_SIM)
from nose.tools import assert_equal as EQ



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

        ARRAY_EQ(c.read(orientation="native", aligned=False),
                 b.read(orientation="native", aligned=False))

        ARRAY_EQ(c.read(orientation="genomic", aligned=False),
                 b.read(orientation="genomic", aligned=False))

        ARRAY_EQ(c.read(orientation="genomic", aligned=True),
                 b.read(orientation="genomic", aligned=True))

        ARRAY_EQ(c.read(orientation="native", aligned=True),
                 b.read(orientation="native", aligned=True))


        # MergeQV is actually not identical presently, due to the
        # clipping at 93 used for ASCII encoding in BAM
        for featureName in [ "DeletionQV",
                             "InsertionQV",
                             "SubstitutionQV",
                             "DeletionTag" ]:
            ARRAY_EQ(c.pulseFeature(featureName, orientation="native", aligned=False),
                     b.pulseFeature(featureName, orientation="native", aligned=False))

        ARRAY_EQ(c.referencePositions(orientation="genomic"),
                 b.referencePositions(orientation="genomic"))

        ARRAY_EQ(c.referencePositions(orientation="native"),
                 b.referencePositions(orientation="native"))

        ARRAY_EQ(c.readPositions(orientation="genomic"),
                 b.readPositions(orientation="genomic"))

        ARRAY_EQ(c.readPositions(orientation="native"),
                 b.readPositions(orientation="native"))


    def testForwardStrandAln(self):
        self.compareAlns(self.cFwd, self.bFwd)


    def testReverseStrandAln(self):
        self.compareAlns(self.cRev, self.bRev)


    def testClipping(self):
        pass


# BT = BamTests()
# BT.test_forward_strand_aln()
# BT.test_reverse_strand_aln()
