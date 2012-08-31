from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
import operator, numpy as np

from GenomicConsensus.plurality.plurality import PluralityWorkerProcess
from AlignmentHitStubs import *

plurality = PluralityWorkerProcess.plurality

def pluralitySequence(refWindow, hits):
    pluralityResults = plurality(refWindow, hits)
    pluralityResults.sort(key=lambda result: result[0])
    locusSummaries = [result[1] for result in pluralityResults]
    return "".join(s.consensus for s in locusSummaries)

def validateStubs(klass):
    actualPluralitySequence = pluralitySequence(klass.referenceWindow,
                                                klass.hits)
    assert_equal(actualPluralitySequence,
                 klass.expectedPluralityConsensus)

class TestPlurality:

    def test_forwardReads(self):
        validateStubs(AllForwardStrandReads)

    def test_reverseReads(self):
        validateStubs(AllReverseStrandReads)

    def test_forwardAndReverseReads(self):
        validateStubs(ForwardAndReverseReads)

    def test_staggeredReads(self):
        validateStubs(StaggeredReads)



