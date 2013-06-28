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
    config = PluralityConfig(minConfidence=0, minCoverage=0)
    css, variants = pluralityConsensusAndVariants(StaggeredReads.referenceWindow,
                                                  StaggeredReads.reference,
                                                  StaggeredReads.hits,
                                                  config)
    assert_equal(StaggeredReads.expectedPluralityConsensus,
                 css.sequence)

    assert_equal(StaggeredReads.expectedPluralityVariants,
                 variants)


def test_computeHaploidVariants():
    config = PluralityConfig(minConfidence=0,
                             minCoverage=0)

    variants1 = _computeVariants(config,
                                 (1, 0, 7),
                                 "GATTACA",
                                 [4]*7,
                                 "GATGACA",
                                 [3]*7,
                                 [35]*7)
    assert_equal([ Variant(1, 3, 4, "T", "G",
                           coverage=4, confidence=35, frequency1=3) ],
                 variants1)

    variants2 = _computeVariants(config,
                                 (1, 0, 7),
                                 "GATTACA",
                                 [4]*7,
                                 ["G", "A", "", "T", "A", "C", "A"],
                                 [3]*7,
                                 [35]*7)
    assert_equal([ Variant(1, 2, 3, "T", "",
                           coverage=4, confidence=35, frequency1=3)],
                 variants2)

    variants3 = _computeVariants(config,
                                 (1, 0, 7),
                                 "GATTACA",
                                 [4]*7,
                                 ["G", "A", "TT", "T", "A", "C", "A"],
                                 [3]*7,
                                 [35]*7)
    assert_equal([ Variant(1, 2, 2, "", "T",
                           coverage=4, confidence=35, frequency1=3)],
                 variants3)


# def test_computeVariantsDiploid():
#     config = PluralityConfig(minConfidence=0,
#                              minCoverage=0,
#                              diploid=True)
#     variants1 = _computeVariants(config,
#                                  (1, 0, 7),
#                                  "GATTACA",
#                                  [20]*7,
#                                  "GATTACA",
#                                  [10]*7,
#                                  [35]*7,
#                                  "GATCACA",
#                                  [10]*7,
#                                  [35]*7)
#     assert_equal([ Substitution(1, 3, 4, "T", "G", 4, 35, 3) ], variants1)
