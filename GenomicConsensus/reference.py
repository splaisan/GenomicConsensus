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

import math, md5, logging, re, numpy as np
from collections import namedtuple, OrderedDict
from pbcore.io import FastaReader, FastaEntry

class ReferenceContig(object):
    """
    A contig from a reference (i.e. FASTA) file.
    """
    def __init__(self, id, name, header, md5sum, sequence, length):
        self.id     = id          # CmpH5-local id
        self.name   = name        # "ref00.."
        self.header = header      # Fasta header
        self.md5sum = md5sum
        self.sequence = sequence
        self.length = length

# Lookup tables.
byName   = OrderedDict()   # Global identifier (string, e.g. "ref000001") -> ReferenceContig
byHeader = OrderedDict()   # Fasta header (string e.g. "chr1") -> ReferenceContig
byId     = OrderedDict()   # CmpH5 local id (integer)    -> ReferenceContig
byMD5    = OrderedDict()   # MD5 sum (e.g. "a13...")     -> ReferenceContig

def idToName(_id):
    return byId[_id].name

def nameToId(name):
    return byName[name].id

def idToHeader(_id):
    return byId[_id].header

def headerToId(header):
    return byHeader[header].id

# Interpret a string key (one of name, header, or id (as string))
# and find the associated id.  Only to be used in interpretation of
# command-line input!
def anyKeyToId(stringKey):
    assert isLoaded()
    if stringKey in byHeader:
        return byHeader[stringKey].id
    elif stringKey in byName:
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
    numFastaEntries = 0
    fastaChecksums = set()
    for fastaEntry in f:
        numFastaEntries += 1
        md5sum = md5.md5(fastaEntry.sequence).hexdigest()
        fastaChecksums.add(md5sum)
        normalizedContigSequence = fastaEntry.sequence.upper()
        if md5sum in cmpH5.referenceTable.MD5:
            cmpH5RefEntry = cmpH5.referenceInfo(md5sum)
            refId         = cmpH5RefEntry.ID
            refName       = cmpH5RefEntry.Name
            contig = ReferenceContig(refId, refName, fastaEntry.raw_name, md5sum,
                                     np.array(normalizedContigSequence, dtype="c"),
                                     len(normalizedContigSequence))
            byId[refId]          = contig
            byName[refName]      = contig
            byMD5[contig.md5sum] = contig
            byHeader[fastaEntry.raw_name] = contig
    logging.info("Loaded %d of %d reference groups from %s " %
                 (len(byId), numFastaEntries, filename))

    # If the cmpH5 has alignments to contigs that weren't contained in
    # the fasta file, report an error.
    cmpH5Checksums = set(cmpH5.referenceTable.MD5)
    if not cmpH5Checksums.issubset(fastaChecksums):
        logging.error("CmpH5 aligned to a contig not represented in FASTA file")
        return 1
    assert isLoaded()

def windowFromString(s):
    assert isLoaded()
    if s == None:
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

def enumerateChunks(refId, referenceStride, referenceWindow=None, overlap=0):
    """
    Enumerate all work chunks (restricted to the window, if provided).
    """
    assert isLoaded()
    assert (referenceWindow == None) or (refId == referenceWindow[0])
    referenceEntry = byId[refId]
    if referenceWindow:
        _, start, end = referenceWindow
    else:
        start, end = (0, referenceEntry.length)

    # The last chunk only needs to reach 'end-overlap', because it will be extended to 'end' -- this prevents the generation of multiple chunks
    # covering the last few bases of the reference (fixes bug #21940)
    for chunkBegin in xrange(start, end-overlap, referenceStride):
        yield (refId,
               max(chunkBegin - overlap, 0),
               min(chunkBegin + referenceStride + overlap, referenceEntry.length))

def numChunks(refId, referenceStride, referenceWindow=None):
    """
    How many chunks will there be for the given refId and window restriction?
    """
    assert isLoaded()
    assert (referenceWindow == None) or (refId == referenceWindow[0])
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
    if referenceWindow == None:
        for refId in byId: yield refId
    else:
        refId, refStart, refEnd = referenceWindow
        yield refId
