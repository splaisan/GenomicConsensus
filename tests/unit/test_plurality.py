from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
import operator, numpy as np

from GenomicConsensus.plurality.plurality import (PluralityConfig,
                                                  pluralityConsensusAndVariants,
                                                  _computeVariants)
from AlignmentHitStubs import *

def test_plurality1():
    css, variants = pluralityConsensusAndVariants(ForwardAndReverseReads.referenceWindow,
                                                  ForwardAndReverseReads.reference,
                                                  ForwardAndReverseReads.hits,
                                                  PluralityConfig())

    assert_equal(ForwardAndReverseReads.expectedPluralityConsensus,
                 css.sequence)

    assert_equal(ForwardAndReverseReads.expectedPluralityVariants,
                 variants)


def test_plurality2():
    css, variants = pluralityConsensusAndVariants(StaggeredReads.referenceWindow,
                                                  StaggeredReads.reference,
                                                  StaggeredReads.hits,
                                                  PluralityConfig(minConfidence=0,
                                                                  minCoverage=0))


    assert_equal(StaggeredReads.expectedPluralityConsensus,
                 css.sequence)

    assert_equal(StaggeredReads.expectedPluralityVariants,
                 variants)


def test_computeVariants():
    variants1 = _computeVariants((1, 0, 7), "GATTACA", "GATGACA", [35]*7, [4]*7, [3]*7)
    assert_equal([ Substitution(1, 3, 4, "T", "G", 4, 35, 3) ], variants1)
    variants2 = _computeVariants((1, 0, 7), "GATTACA",
                                 ["G", "A", "", "T", "A", "C", "A"],
                                 [35]*7, [4]*7, [3]*7)
    assert_equal([Deletion(1, 2, 3, "T", "", 4, 35, 3)], variants2)
    variants2 = _computeVariants((1, 0, 7), "GATTACA",
                                 ["G", "A", "TT", "T", "A", "C", "A"],
                                 [35]*7, [4]*7, [3]*7)
    assert_equal([Insertion(1, 2, 2, "", "T", 4, 35, 3)], variants2)
