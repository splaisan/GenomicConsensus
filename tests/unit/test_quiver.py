from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal

from GenomicConsensus.quiver import utils
from GenomicConsensus.variants import *
from ConsensusCore import *

class TestVariantsFromAlignment(object):

    def testVariantsFromAlignment1(self):
        a = PairwiseAlignment("GATTACA", "GAT-ACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Deletion(1, 1003, 1004, "T", "") ], vs)

    def testVariantsFromAlignment2(self):
        a = PairwiseAlignment("GA-TTACA", "GATTTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Insertion(1, 1002, 1002, "", "T") ], vs)

    def testVariantsFromAlignment3(self):
        a = PairwiseAlignment("GATTACA", "GAGGACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Substitution(1, 1002, 1004, "TT", "GG") ], vs)

    def testNoCallBasesInReference1(self):
        a = PairwiseAlignment("GATTNGATT", "GAGGATATT")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Substitution(1, 1002, 1004, "TT", "GG"),
                       Substitution(1, 1005, 1006, "G",  "T")   ], vs)

    def testLongInsertion(self):
        a = PairwiseAlignment("GA-----TTACA", "GACCCCCTTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Insertion(1, 1002, 1002, "", "CCCCC") ], vs)

    def testLongDeletion(self):
        a = PairwiseAlignment("GACCCCCTTACA", "GA-----TTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Deletion(1, 1002, 1007, "CCCCC", "") ], vs)

    def testTwoSubstitutions(self):
        a = PairwiseAlignment("GATTACA", "GAGTAGA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Substitution(1, 1002, 1003, "T", "G"),
                       Substitution(1, 1005, 1006, "C", "G") ], vs)
