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

    # Is this the behavior we want?  We could output multibase
    # substitutions, too.
    def testVariantsFromAlignment3(self):
        a = PairwiseAlignment("GATTACA", "GAGGACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Substitution(1, 1002, 1004, "TT", "GG") ], vs)

    def testLongInsertion(self):
        a = PairwiseAlignment("GA-----TTACA", "GACCCCCTTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Insertion(1, 1002, 1002, "", "CCCCC") ], vs)

    def testLongDeletion(self):
        a = PairwiseAlignment("GACCCCCTTACA", "GA-----TTACA")
        vs = utils.variantsFromAlignment(a, (1, 1000, 2000))
        assert_equal([ Deletion(1, 1002, 1007, "CCCCC", "") ], vs)


class TestInverseMutations(object):

    def testMutationEquality(self):
        # Make sure Mutation equality operator works correctly
        assert_equal(Mutation(SUBSTITUTION, 10, "T"),
                     Mutation(SUBSTITUTION, 10, "T"))

    def testInverseOfSubstitution(self):
        inv = utils.inverseMutations(1000, Substitution(1, 1010, 1011, "T", "G"))
        assert_equal([Mutation(SUBSTITUTION, 10, "T")], inv)

    def testInverseOfInsertion(self):
        inv = utils.inverseMutations(1000, Insertion(1, 1010, 1010, "", "G"))
        assert_equal([Mutation(DELETION, 10, "-")], inv)

    def testInverseOfDeletion(self):
        inv = utils.inverseMutations(1000, Deletion(1, 1010, 1011, "G", ""))
        assert_equal([Mutation(INSERTION, 10, "G")], inv)

    def testInverseOfLongInsertion(self):
        inv = utils.inverseMutations(1000, Insertion(1, 1010, 1010, "", "GAT"))
        assert_equal([Mutation(DELETION, 10, "-"),
                      Mutation(DELETION, 11, "-"),
                      Mutation(DELETION, 12, "-")], inv)

    def testInverseOfLongDeletion(self):
        inv = utils.inverseMutations(1000, Deletion(1, 1010, 1013, "GAT", ""))
        assert_equal([Mutation(INSERTION, 10, "G"),
                      Mutation(INSERTION, 10, "A"),
                      Mutation(INSERTION, 10, "T")], inv)

    def testInverseOfLongSubstitution(self):
        inv = utils.inverseMutations(1000, Substitution(1, 1010, 1013, "GGG", "CCC"))
        assert_equal([Mutation(SUBSTITUTION, 10, "G"),
                      Mutation(SUBSTITUTION, 11, "G"),
                      Mutation(SUBSTITUTION, 12, "G")], inv)

