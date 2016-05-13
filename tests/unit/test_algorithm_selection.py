
from nose.tools import assert_equal as EQ
from GenomicConsensus.algorithmSelection import bestAlgorithm_


def test_algorithm_selection():
    EQ("quiver", bestAlgorithm_(["P6-C4"]))
    EQ("quiver", bestAlgorithm_(["P6-C4", "P5-C3"]))
    EQ("arrow",  bestAlgorithm_(["S/P1-C1/beta"]))
    EQ("arrow",  bestAlgorithm_(["P6-C4", "S/P1-C1/beta"]))
    EQ(None,     bestAlgorithm_(["P6-C4", "unknown"]))
    EQ("quiver", bestAlgorithm_(["S/P1-C1"]))
    EQ("arrow",  bestAlgorithm_(["P6-C4", "S/P1-C1"]))
    EQ("arrow",  bestAlgorithm_(["P5-C3", "S/P1-C1"])) # (Arrow pres. no training for P5.  But it will tell us that)
