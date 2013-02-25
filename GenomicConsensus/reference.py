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
from pbcore.io import FastaReader, FastaRecord

class ReferenceContig(object):
    """
    A contig from a reference (i.e. FASTA) file.
    """
    def __init__(self, id, name, md5sum, sequence, length):
        self.id       = id          # CmpH5-local id
        self.name     = name        # Fasta header
        self.md5sum   = md5sum
        self.sequence = sequence
        self.length   = length

byName   = OrderedDict()   # Fasta header (string e.g. "chr1") -> FastaRecord
byId     = OrderedDict()   # CmpH5 local id (integer)          -> FastaRecord
byMD5    = OrderedDict()   # MD5 sum (e.g. "a13...")           -> FastaRecord

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

def isLoaded():
    return bool(byMD5)

def loadFromFile(filename, cmpH5):
    """
    Reads reference from FASTA file, loading
    lookup tables that can be used any time later.
    """
    # Load contigs
    assert not isLoaded()
    f = FastaReader(filename)
    numFastaRecords = 0
    fastaChecksums = set()
    for fastaRecord in f:
        numFastaRecords += 1
        md5sum = fastaRecord.md5
        fastaChecksums.add(md5sum)
        normalizedContigSequence = fastaRecord.sequence.upper()
        if md5sum in cmpH5.referenceInfoTable.MD5:
            cmpH5RefEntry = cmpH5.referenceInfo(md5sum)
            refId         = cmpH5RefEntry.ID
            refName       = fastaRecord.name
            contig = ReferenceContig(refId, refName, md5sum,
                                     np.array(normalizedContigSequence, dtype="c"),
                                     len(normalizedContigSequence))
            byId[refId]          = contig
            byName[refName]      = contig
            byMD5[contig.md5sum] = contig
    logging.info("Loaded %d of %d reference groups from %s " %
                 (len(byId), numFastaRecords, filename))

    # If the cmpH5 has alignments to contigs that weren't contained in
    # the fasta file, report an error.
    cmpH5Checksums = set(cmpH5.referenceInfoTable.MD5)
    if not cmpH5Checksums.issubset(fastaChecksums):
        logging.error("CmpH5 aligned to a contig not represented in FASTA file")
        return 1
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


def enumerateChunks(refId, referenceStride, referenceWindow=None, overlap=0):
    """
    Enumerate all work chunks (restricted to the window, if provided).
    """
    assert isLoaded()
    assert (referenceWindow is None) or (refId == referenceWindow[0])
    referenceEntry = byId[refId]
    if referenceWindow:
        _, start, end = referenceWindow
    else:
        start, end = (0, referenceEntry.length)

    # The last chunk only needs to reach 'end-overlap', because it
    # will be extended to 'end' -- this prevents the generation of
    # multiple chunks covering the last few bases of the reference
    # (fixes bug #21940)
    for chunkBegin in xrange(start, end-overlap, referenceStride):
        yield (refId,
               max(chunkBegin - overlap, 0),
               min(chunkBegin + referenceStride + overlap, referenceEntry.length))

def numChunks(refId, referenceStride, referenceWindow=None):
    """
    How many chunks will there be for the given refId and window restriction?
    """
    assert isLoaded()
    assert (referenceWindow is None) or (refId == referenceWindow[0])
    referenceEntry = byId[refId]
    if referenceWindow:
        _, start, end = referenceWindow
    else:
        start, end = (0, referenceEntry.length)
    return int(math.ceil(float(end-start)/referenceStride))

def enumerateIds(referenceWindow=None):
    """
    Enumerate all refIds (subject to the referenceWindow restriction, if provided).
    """
    assert isLoaded()
    if referenceWindow is None:
        for refId in byId: yield refId
    else:
        refId, refStart, refEnd = referenceWindow
        yield refId
