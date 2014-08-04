
from nose.tools import assert_equal as EQ
from GenomicConsensus.quiver.utils import findDinucleotideRepeats


def test_findDinucleotideRepeats():
    EQ([], findDinucleotideRepeats("GATTACA"))
    EQ([((1,7), "AT")], findDinucleotideRepeats("GATATATACA"))
    EQ([((1, 7), "AT"), ((7, 13), "AC")],
       findDinucleotideRepeats("GATATATACACACA"))
