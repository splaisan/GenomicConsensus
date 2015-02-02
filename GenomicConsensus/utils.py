#################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

# Author: David Alexander

from __future__ import absolute_import
import math, numpy as np, os.path, sys, itertools

def die(msg):
    print msg
    sys.exit(-1)

class CommonEqualityMixin(object):
    def __eq__(self, other):
        return (isinstance(other, self.__class__)
            and self.__dict__ == other.__dict__)

    def __ne__(self, other):
        return not self.__eq__(other)


# An exception for incompatible cmp.h5 files
class IncompatibleDataException(Exception):
    pass

# We truncate QVs at 93 because the FASTQ format downstream can only
# support QVs in the range [0, 93] without lossage.

def error_probability_to_qv(error_probability, cap=93):
    """
    Convert an error probability to a phred-scaled QV.
    """
    if error_probability==0:
        return cap
    else:
        return min(cap, int(round(-10*math.log10(error_probability))))


_complement = { "A" : "T",
                "C" : "G",
                "G" : "C",
                "T" : "A",
                "-" : "-" }

def complement(s):
    cStr = "".join(_complement[c] for c in s)
    if type(s) == str:
        return cStr
    else:
        return np.fromstring(cStr, "S1")

def reverseComplement(s):
    return complement(s)[::-1]

def fileFormat(filename):
    if filename.endswith(".gz"):
        ext = os.path.splitext(filename[:-3])[1]
    else:
        ext = os.path.splitext(filename)[1]
    ext = ext.lower()
    if   ext in [".fa", ".fasta"]: return "FASTA"
    elif ext in [".fq", ".fastq"]: return "FASTQ"
    elif ext in [".gff" ]:         return "GFF"
    elif ext in [".csv" ]:         return "CSV"
    else: raise Exception, "Unrecognized file format"

def rowNumberIsInReadStratum(readStratum, rowNumber):
    n, N = readStratum
    return (rowNumber % N) == n

def readsInWindow(cmpH5, window, depthLimit=None,
                  minMapQV=0, strategy="fileorder",
                  stratum=None, barcode=None):
    """
    Return up to `depthLimit` reads (as row numbers integers) where
    the mapped reference intersects the window.  If depthLimit is None,
    return all the reads meeting the criteria.

    `strategy` can be:
      - "longest" --- get the reads with the longest length in the window
      - "spanning" --- get only the reads spanning the window
      - "fileorder" --- get the reads in file order
    """
    assert strategy in {"longest", "spanning", "fileorder"}

    if stratum is not None:
        raise ValueError, "stratum needs to be reimplemented"

    def depthCap(iter):
        if depthLimit is not None:
            return list(itertools.islice(iter, 0, depthLimit))
        else:
            return list(iter)

    def lengthInWindow(hit):
        return min(hit.tEnd, winEnd) - max(hit.tStart, winStart)

    winId, winStart, winEnd = window
    if barcode == None:
        alnHits = ( hit
                    for hit in cmpH5.readsInRange(winId, winStart, winEnd)
                    if hit.MapQV >= minMapQV )
    else:
        alnHits = ( hit
                    for hit in cmpH5.readsInRange(winId, winStart, winEnd)
                    if ((hit.MapQV >= minMapQV) and
                        (hit.barcode == barcode)) )

    if strategy == "fileorder":
        return depthCap(alnHits)
    elif strategy == "spanning":
        winLen = winEnd - winStart
        return depthCap( hit for hit in alnHits
                         if lengthInWindow(hit) == winLen )
    elif strategy == "longest":
        # Well, this defeats the iterable/laziness, performs poorly if
        # deep coverage.  sorted is stable.
        return depthCap(sorted(alnHits, key=lengthInWindow, reverse=True))


def datasetCountExceedsThreshold(cmpH5, threshold):
    """
    Does the file contain more than `threshold` datasets?  This
    impacts whether or not we should disable the chunk cache.
    """
    total = 0
    for i in np.unique(cmpH5.AlnGroupID):
        total += len(cmpH5._alignmentGroup(i))
        if total > threshold:
            return True
    return False

#
# Some lisp functions we want
#
fst   = lambda t: t[0]
snd   = lambda t: t[1]
third = lambda t: t[2]

def nub(it):
    """
    Unique entries in an iterable, preserving order
    """
    seen = set()
    for x in it:
        if x not in seen: yield(x)
        seen.add(x)
