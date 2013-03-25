from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal

from GenomicConsensus.consensus import *

def test_areContiguous():
    refWindowsGood = [ (0, 100, 200),
                       (0, 200, 300) ]
    refWindowsBad1 = [ (0, 100, 200),
                       (0, 201, 300) ]
    refWindowsBad2 = [ (0, 100, 200),
                       (1, 200, 300) ]

    assert areContiguous(refWindowsGood)
    assert not areContiguous(refWindowsBad1)

def test_join():
    chunks = [ Consensus( (0, 100, 110), "GATTACA", range(7) ),
               Consensus( (0, 110, 120), "CATTACA", range(7,14) ) ]
    joined = join(chunks)
    expectedJoined = Consensus( (0, 100, 120), "GATTACACATTACA", range(14))
    assert_equal(expectedJoined, joined)
