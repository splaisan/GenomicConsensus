
from nose.tools import assert_equal as EQ
from GenomicConsensus.algorithmSelection import bestAlgorithm


def test_algorithm_selection():
    EQ("quiver", bestAlgorithm(["P6-C4"]))
    EQ("quiver", bestAlgorithm(["P6-C4", "P5-C3"]))
    EQ("poa",    bestAlgorithm(["S/P1-C1/beta"]))
    EQ("poa",    bestAlgorithm(["P6-C4", "S/P1-C1/beta"]))
    EQ("arrow",  bestAlgorithm(["S/P1-C1"]))
    EQ("arrow",  bestAlgorithm(["P6-C4", "S/P1-C1"]))
    EQ("arrow",  bestAlgorithm(["P5-C3", "S/P1-C1"])) # (Arrow pres. no training for P5.  But it will tell us that)
