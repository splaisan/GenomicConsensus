from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal

from GenomicConsensus.io.fastx import FastaWriter, FastqWriter
from StringIO import StringIO
import numpy as np

def testFastaWriter():
    s = StringIO()
    w = FastaWriter(s)

    w.writeRecord("foo", "GATTACA"*20)
    w.writeRecord("bar", "")
    w.writeRecord("baz", "T"*60)

    expected = """\
>foo
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT
ACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG
ATTACAGATTACAGATTACA
>bar
>baz
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"""
    assert_equal(expected, s.getvalue())
    w.close()

def qvArray(arr):
    return np.array(arr, dtype=np.uint8)

def testFastqWriter():
    s = StringIO()
    w = FastqWriter(s)

    w.writeRecord("foo", "GATTACA"*20, qvArray([1,2,3,4,5,6,7]*20))
    w.writeRecord("bar", "", qvArray([]))
    w.writeRecord("baz", "T"*60, qvArray(range(60, 120)))


    expected = """\
@foo
GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATT
ACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAG
ATTACAGATTACAGATTACA
+
"#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%
&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("#$%&'("
#$%&'("#$%&'("#$%&'(
@bar
+
@baz
TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
+
]^_`abcdefghijklmnopqrstuvwxyz{|}~~~~~~~~~~~~~~~~~~~~~~~~~~~"""

    print expected
    assert_equal(expected, s.getvalue())
    w.close()
