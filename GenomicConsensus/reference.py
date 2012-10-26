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
byName  = OrderedDict()   # Global identifier (string, e.g. "ref000001") -> ReferenceContig
byId    = OrderedDict()   # CmpH5 local id (integer)    -> ReferenceContig
byMD5   = OrderedDict()   # MD5 sum (e.g. "a13...")     -> ReferenceContig

def idToName(_id):
    return byId[_id].name

def nameToId(name):
    return byName[name].id

def idToHeader(_id):
    return byId[_id].header

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
        refId    = nameToId(m.group(1))
        refStart = int(m.group(2))
        refEnd   = min(int(m.group(3)), byId[refId].length)
    else:
        refId    = nameToId(s)
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
    for chunkBegin in xrange(start, end, referenceStride):
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
