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

import math, logging, re, numpy as np
from collections import OrderedDict
from pbcore.io import FastaTable

from ConsensusCore import CoveredIntervals
from .utils import holes, nub, die
from .Chunk import WorkChunk


class UppercasingMmappedFastaSequence(object):

    def __init__(self, mmappedFastaSequence):
        self.other = mmappedFastaSequence

    def __getitem__(self, spec):
        snip = self.other.__getitem__(spec)
        return snip.upper()

class ReferenceContig(object):
    """
    A contig from a reference (i.e. FASTA) file.
    """
    def __init__(self, id, name, sequence, length):
        self.id        = id          # CmpH5-local id
        self.name      = name        # Fasta header
        self.sequence  = UppercasingMmappedFastaSequence(sequence)
        self.length    = length


byName   = OrderedDict()   # Fasta header (string e.g. "chr1") -> FastaRecord
byId     = OrderedDict()   # CmpH5 local id (integer)          -> FastaRecord

def idToName(_id):
    return byId[_id].name

def nameToId(name):
    return byName[name].id

# Interpret a string key (one of name, or id (as string))
# and find the associated id.  Only to be used in interpretation of
# command-line input!
def anyKeyToId(stringKey):
    assert isLoaded()
    if stringKey in byName:
        return byName[stringKey].id
    elif stringKey.isdigit():
        refId = int(stringKey)
        return byId[refId].id
    else:
        raise Exception, "Unknown reference name: %s" % stringKey

def sequenceInWindow(window):
    refId, refStart, refEnd = window
    return byId[refId].sequence[refStart:refEnd]

filename = None

def isLoaded():
    return filename != None

def loadFromFile(filename_, cmpH5):
    """
    Reads reference from FASTA file, loading
    lookup tables that can be used any time later.
    """
    # Contigs in FASTA may disagree with those in cmp.h5 ref info
    # table, for instance if the FASTA has been edited.  Here's how we
    # handle things:
    #
    # |fastaContigs \   cmpContigs| > 0 : OK, extra FASTA contigs just ignored
    # |cmpContigs   \ fastaContigs| > 0 : Not necessarily OK---a warning should be
    #                                     issued.  We then proceed to operate on
    #                                     the contigs that are in both.
    # |cmpContigs ^ fastaContigs| == 0  : Nothing to work with.  This is an error.
    #
    # While we formerly used MD5s to vouch for the identity of a
    # contig, we not use the name.  This is an inferior approach but
    # is necessary, in using the FastaTable.

    # Load contigs
    assert not isLoaded()
    try:
        f = FastaTable(filename_)
    except IOError as e:
        die(e)

    cmpContigNames = set(cmpH5.referenceInfoTable.FullName)

    for fastaRecord in f:
        refName = fastaRecord.name
        if refName in cmpContigNames:
            cmpH5RefEntry = cmpH5.referenceInfo(refName)
            refId         = cmpH5RefEntry.ID
            sequence      = UppercasingMmappedFastaSequence(fastaRecord.sequence)
            length        = len(fastaRecord.sequence)
            contig = ReferenceContig(refId, refName, sequence, length)
            byId[refId]          = contig
            byName[refName]      = contig
    loadedFastaContigNames = set(byName.keys())
    logging.info("Loaded %d of %d reference groups from %s " %
                 (len(byId), len(loadedFastaContigNames), filename_))

    if len(byId) == 0:
        die("No reference groups in the FASTA file were aligned against.  " \
            "Did you select the wrong reference FASTA file?")
    elif (cmpContigNames - loadedFastaContigNames):
        logging.warn(
            "Some reference contigs aligned against are not found in " \
            "the reference FASTA.  Will process only those contigs "   \
            "supported by the reference FASTA.")

    global filename
    filename = filename_
    assert isLoaded()

def stringToWindow(s):
    assert isLoaded()
    if s is None:
        return None
    m = re.match("(.*):(.*)-(.*)", s)
    if m:
        refId    = anyKeyToId(m.group(1))
        refStart = int(m.group(2))
        refEnd   = min(int(m.group(3)), byId[refId].length)
    else:
        refId    = anyKeyToId(s)
        refStart = 0
        refEnd   = byId[refId].length
    return (refId, refStart, refEnd)

def windowToString(referenceWindow):
    assert isLoaded()
    refId, refStart, refEnd = referenceWindow
    return "%s:%d-%d" % (idToName(refId),
                         refStart,
                         refEnd)

def enumerateSpans(refId, referenceWindows=()):
    """
    Enumerate the contiguous spans along this reference contig that
    are to be analyzed.
    """
    assert isLoaded()
    referenceEntry = byId[refId]
    referenceEntrySpan = (refId, 0, referenceEntry.length)

    for refWin in (referenceWindows or [referenceEntrySpan]):
        refWinId, start, end = refWin
        if refWinId == refId:
            yield (refId, start, end)

def enumerateChunks(refId, referenceStride, referenceWindows=()):
    """
    Enumerate all work chunks on this reference contig (restricted to
    the windows, if provided).
    """
    for span in enumerateSpans(refId, referenceWindows):
        _, start, end = span
        for chunkBegin in xrange(start, end, referenceStride):
            win = (refId,
                   chunkBegin,
                   min(chunkBegin + referenceStride, end))
            yield WorkChunk(win, True)

def fancyEnumerateChunks(cmpH5, refId, referenceStride,
                         minCoverage, minMapQV, referenceWindows=()):
    """
    Enumerate chunks, creating chunks with hasCoverage=False for
    coverage cutouts.
    """
    startRow = cmpH5.referenceInfo(refId).StartRow
    endRow   = cmpH5.referenceInfo(refId).EndRow
    goodMapQVs = (cmpH5.MapQV[startRow:endRow] >= minMapQV)
    tStart = cmpH5.tStart[startRow:endRow][goodMapQVs]
    tEnd   = cmpH5.tEnd  [startRow:endRow][goodMapQVs]

    for span in enumerateSpans(refId, referenceWindows):
        _, spanStart, spanEnd = span
        coveredIntervals = kCoveredIntervals(minCoverage, tStart, tEnd, spanStart, spanEnd)
        unCoveredIntervals = holes(span, coveredIntervals)

        for (s, e) in sorted(list(coveredIntervals) + unCoveredIntervals):
            win = (refId, s, e)
            if (s, e) in coveredIntervals:
                for chunk in enumerateChunks(refId, referenceStride, [(refId, s, e)]):
                    yield chunk
            else:
                yield WorkChunk(win, False)


def numReferenceBases(refId, referenceWindows=()):
    """
    Termination is determined to be when the result collector has
    built consensus corresponding to the exact number of reference
    bases in the window under consideration.
    """
    return sum((end - start)
               for (_, start, end) in enumerateSpans(refId, referenceWindows))


def enumerateIds(referenceWindows=()):
    """
    Enumerate all refIds (subject to the referenceWindows restriction,
    if provided).
    """
    assert isLoaded()
    if referenceWindows == ():
        for refId in byId:
            yield refId
    else:
        for refId in nub(refId for (refId, _, _) in referenceWindows):
            yield refId

def enlargedReferenceWindow(refWin, overlap):
    assert isLoaded()
    refId, refStart, refEnd = refWin
    contigLength = byId[refId].length
    return (refId,
            max(0, refStart - overlap),
            min(refEnd + overlap, contigLength))

def kCoveredIntervals(k, tStart, tEnd, winStart, winEnd):
    return CoveredIntervals(k, tStart, tEnd, int(winStart), int(winEnd-winStart))
