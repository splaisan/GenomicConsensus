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
import math, numpy as np, os.path, sys
from pbcore.io.rangeQueries import projectIntoRange

def die(msg):
    print msg
    sys.exit(-1)

# Some lisp functions we want
fst   = lambda t: t[0]
snd   = lambda t: t[1]
third = lambda t: t[2]

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

def probability_to_qv(probability, cap=93):
    """
    Convert an event probability (not an error probability) to a
    phred-scaled QV.
    """
    return error_probability_to_qv(1.0-probability, cap)


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

def noEvidenceConsensusCall(referenceSequence, noEvidenceConsensusCallMode):
    # referenceSequence is str, return value is str
    # return value type is the same type as input
    assert isinstance(referenceSequence, str)

    if noEvidenceConsensusCallMode == "reference":
        result = referenceSequence
    elif noEvidenceConsensusCallMode == "lowercasereference":
        result =  referenceSequence.lower()
    elif noEvidenceConsensusCallMode == "nocall":
        result = "N" * len(referenceSequence)
    else:
        raise Exception, "Invalid `noEvidenceConsensusCallMode`"

    return result

class Window(object):
    def __init__(self, contigKey, start, end):
        self.contigKey  = contigKey
        self.start = start
        self.end   = end

    def __repr__(self):
        return "%s:%d-%d" % (self.contigKey, self.start, self.end)

    def __len__(self):
        return (self.end - self.start)

    def __eq__(self, other):
        return (self.contigKey == other.contigKey and
                self.start     == other.start and
                self.end       == other.end)

    @staticmethod
    def fromString(windowString):
        splitOnColon = windowString.split(":")
        if len(splitOnColon) < 2:
            key = 1
        else:
            key, rest = splitOnColon
        if "-" not in rest:
            refPos = int(rest)
            start = refPos - 10
            end   = refPos + 10
        else:
            start, end = map(int, rest.split("-"))
        assert end >= start
        return Window(key, start, end)

    @staticmethod
    def fromGffString(windowString):
        w = Window.fromString(windowString)
        w.start -= 1
        return w

    def subWindow(self, start, end):
        assert start >= self.start
        assert end   <= self.end
        return Window(self.contigKey, start, end)


def readsInWindow(cmpH5, window, depthLimit=None, minMapQV=0, strategy="fileorder"):
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

    def depthCap(lst):
        if depthLimit is not None: return lst[:depthLimit]
        else: return lst

    contigLocalId = cmpH5.referenceInfo(window.contigKey).ID

    rowNumbers = cmpH5.readsInRange(contigLocalId, window.start, window.end,
                                    justIndices=True)

    rowNumbers = rowNumbers[cmpH5.MapQV[rowNumbers] >= minMapQV]

    if strategy == "fileorder":
        return depthCap(rowNumbers)

    tStartTruncated = np.maximum(window.start, cmpH5.tStart[rowNumbers])
    tEndTruncated   = np.minimum(window.end,   cmpH5.tEnd[rowNumbers])
    lengthsInWindow = tEndTruncated - tStartTruncated

    if strategy == "spanning":
        return depthCap(rowNumbers[lengthsInWindow==(window.end-window.start)])
    elif strategy == "longest":
        ordering = np.lexsort((rowNumbers, -lengthsInWindow))
        return depthCap(rowNumbers[ordering])


def kSpannedIntervals(refWindow, k, start, end):
    """
    Find intervals in the window that are k-spanned by the reads.

    Given:
     `refWindow`: the window under consideration
     `k`: the number of reads that must span intervals to be returned
     `start`, `end`: numpy arrays of start and end coordinates for reads,
       where the extent of each read is [start, end).  Must be ordered
       so that `start` is sorted in ascending order.

    Find a maximal set of maximal disjoint intervals within
    refWindow such that each interval is spanned by at least k reads.
    Intervals are returned in sorted order, as a list of (start, end)
    tuples.

    Note that this is a greedy search procedure and may not always
    return the optimal solution, in some sense.  However it will
    always return the optimal solutions in the most common cases.
    """
    assert k >= 1

    # Translate the start, end to coordinate system where
    # refWindow.start is 0.
    start = start - refWindow.start
    end   = end - refWindow.start
    winStart = 0
    winEnd   = refWindow.end - refWindow.start

    # Truncate to bounds implied by refWindow
    start = np.maximum(winStart, start)
    end   = np.minimum(winEnd,   end)

    positions = np.arange(winEnd - winStart, dtype=int)
    coverage = projectIntoRange(start, end,
                                winStart, winEnd)
    x = -1
    y = 0
    intervalsFound = []

    while y < winEnd:
        # Step 1: let x be the first pos >= y that is k-covered
        eligible = np.flatnonzero((positions >= y) & (coverage >= k))
        if len(eligible) > 0:
            x = eligible[0]
        else:
            break

        # Step 2: extend the window [x, y) until [x, y) is no longer
        # k-spanned.  Do this by setting y to the k-th largest `end`
        # among reads covering x
        eligible = end[(start <= x)]
        eligible.sort()
        if len(eligible) >= k:
            y = eligible[-k]
        else:
            break

        intervalsFound.append((x, y))

    # Translate intervals back
    return [ (s + refWindow.start,
              e + refWindow.start) for (s, e) in intervalsFound ]


def abut(intervals):
    """
    Abut adjacent intervals.  Useful for debugging...
    """
    output = []
    lastS = None
    lastE = None
    for (s, e) in intervals:
        if s == lastE:
            lastS, lastE = lastS, e
        else:
            if lastS is not None:
                output.append((lastS, lastE))
            lastS, lastE = s, e
    output.append((lastS, lastE))
    return output


def holes(refWindow, intervals):
    """
    Given a window and a set of disjoint subintervals, return the
    "holes", which are the intervals of the refWindow not covered by
    the given subintervals.
    """
    output = []
    intervals = sorted(intervals)
    lastE = refWindow.start
    for (s, e) in intervals:
        if s > lastE:
            output.append((lastE, s))
        lastE = e
    if lastE < refWindow.end:
        output.append((lastE, refWindow.end))
    return output
