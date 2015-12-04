from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal

from GenomicConsensus.quiver import utils
from GenomicConsensus.variants import *
from ConsensusCore import *

class TestVariantsFromAlignment(object):

    def testVariantsFromAlignment1(self):
        a = PairwiseAlignment("GATTACA", "GAT-ACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Variant(1, 1003, 1004, "T", "") ], vs)

    def testVariantsFromAlignment2(self):
        a = PairwiseAlignment("GA-TTACA", "GATTTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Variant(1, 1002, 1002, "", "T") ], vs)

    def testVariantsFromAlignment3(self):
        a = PairwiseAlignment("GATTACA", "GAGGACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Variant(1, 1002, 1004, "TT", "GG") ], vs)

    def testVariantsFromAlignment4(self):
        a = PairwiseAlignment("GA-TACA", "GATTACA")
        qvs = [0, 0, 1, 0, 0, 0, 0]
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000), qvs)
        assert_equal([ Variant(1, 1002, 1002, "", "T", confidence=1) ], vs)

    def testVariantsFromAlignment5(self):
        a = PairwiseAlignment("-ATTACA", "GATTACA")
        qvs = [1, 0, 0, 0, 0, 0, 0]
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000), qvs)
        assert_equal([ Variant(1, 1000, 1000, "", "G", confidence=1) ], vs)

    def testVariantsFromAlignment6(self):
        a = PairwiseAlignment("GATTAC-", "GATTACA")
        qvs = [0, 0, 0, 0, 0, 0, 1]
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000), qvs)
        assert_equal([ Variant(1, 1006, 1006, "", "A", confidence=1) ], vs)

    # this test is completely nonsensical, deletion events have no meaningful
    # confidence represented in the QV track. But this is the expected behavior
    def testVariantsFromAlignment7(self):
        a = PairwiseAlignment("GATTACA", "GATTAC-")
        qvs = [0, 0, 0, 0, 0, 1]
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000), qvs)
        assert_equal([ Variant(1, 1006, 1007, "A", "", confidence=1) ], vs)

    # also totally bogus for similar reasons
    def testVariantsFromAlignment8(self):
        a = PairwiseAlignment("GATTACA", "-ATTACA")
        qvs = [1, 0, 0, 0, 0, 0]
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000), qvs)
        assert_equal([ Variant(1, 1000, 1001, "G", "", confidence=1) ], vs)

    def testNoCallBasesInReference1(self):
        a = PairwiseAlignment("GATTNGATT", "GAGGATATT")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Variant(1, 1002, 1004, "TT", "GG"),
                       Variant(1, 1005, 1006, "G",  "T")   ], vs)

    def testLongInsertion(self):
        a = PairwiseAlignment("GA-----TTACA", "GACCCCCTTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Variant(1, 1002, 1002, "", "CCCCC") ], vs)

    def testLongDeletion(self):
        a = PairwiseAlignment("GACCCCCTTACA", "GA-----TTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Variant(1, 1002, 1007, "CCCCC", "") ], vs)

    def testTwoSubstitutions(self):
        a = PairwiseAlignment("GATTACA", "GAGTAGA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Variant(1, 1002, 1003, "T", "G"),
                       Variant(1, 1005, 1006, "C", "G") ], vs)
