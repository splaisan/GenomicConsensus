
__all__ = ["BamReader", "BamAlignment"]

from pysam import Samfile
from pbcore.io import BasH5Collection
from pbcore.io._utils import rec_join  # FIXME


import numpy as np
from itertools import groupby

class UnavailableFeature(Exception): pass
class Unimplemented(Exception): pass

class BamReader(object):
    def __init__(self, fname):
        self.filename = fname
        self.dual = Samfile(fname, "rb")

        # Make ref table.
        refRecords = self.dual.header["SQ"]
        refNames   = [r["SN"] for r in refRecords]
        refLengths = [r["LN"] for r in refRecords]
        refMD5s    = [r["M5"] for r in refRecords]
        refIds = map(self.dual.gettid, refNames)
        nRefs = len(refRecords)

        self._referenceInfoTable = np.rec.fromrecords(zip(
            refIds,
            refIds,
            refNames,
            refNames,
            refLengths,
            refMD5s,
            np.zeros(nRefs, dtype=np.uint32),
            np.zeros(nRefs, dtype=np.uint32)),
            dtype=[('ID', '<i8'), ('RefInfoID', '<i8'),
                   ('Name', 'O'), ('FullName', 'O'),
                   ('Length', '<i8'), ('MD5', 'O'),
                   ('StartRow', '<u4'), ('EndRow', '<u4')])
        self._referenceDict = {}
        self._referenceDict.update(zip(refIds, self._referenceInfoTable))
        self._referenceDict.update(zip(refNames, self._referenceInfoTable))

        # Make movie table

    def attach(self, fofnFilename):
        self.basH5Collection = BasH5Collection(fofnFilename)

    @property
    def moviesAttached(self):
        return (self.basH5Collection is not None)

    @property
    def alignmentIndex(self):
        raise UnavailableFeature("BAM has no alignment index")

    @property
    def movieInfoTable(self):
        raise Unimplemented()

    @property
    def referenceInfoTable(self):
        return self._referenceInfoTable

    @property
    def readType(self):
        raise Unimplemented()

    @property
    def version(self):
        raise Unimplemented()

    def versionAtLeast(self, minimalVersion):
        raise Unimplemented()

    @property
    def primaryVersion(self):
        raise Unimplemented()

    def primaryVersionAtLeast(self, minimalVersion):
        raise Unimplemented()

    def softwareVersion(self, programName):
        raise Unimplemented()

    @property
    def isSorted(self):
        return True

    @property
    def isBarcoded(self):
        raise Unimplemented()

    @property
    def isEmpty(self):
        return (len(self) == 0)

    def alignmentGroup(self, alnGroupId):
        raise UnavailableFeature("BAM has no HDF5 groups")

    @property
    def movieNames(self):
        raise Unimplemented()

    def movieInfo(self, movieId):
        raise Unimplemented()

    def referenceInfo(self, key):
        return self._referenceDict[key]

    def readsInRange(self, winId, winStart, winEnd, justIndices=False):
        # PYSAM BUG: fetch doesn't work if arg 1 is tid and not rname
        if not isinstance(winId, str):
            winId = self.dual.getrname(winId)
        if justIndices == True:
            raise UnavailableFeature("BAM is not random-access")
        else:
            return ( BamAlignment(self, it)
                     for it in self.dual.fetch(winId, winStart, winEnd) )

    def hasPulseFeature(self, featureName):
        return False

    def pulseFeaturesAvailable(self):
        return []

    @property
    def barcode(self):
        raise Unimplemented()

    @property
    def barcodeName(self):
        raise Unimplemented()

    @property
    def barcodes(self):
        raise Unimplemented()

    def __repr__(self):
        return "<BamReader for %s>" % self.filename

    def __getitem__(self, rowNumbers):
        raise UnavailableFeature("BAM doesn't support true random access")

    def __iter__(self):
        for a in self.dual:
            yield BamAlignment(self, a)

    def __len__(self):
        return self.dual.mapped

    # def __getattr__(self, key):
    #     raise Unimplemented()

    def close(self):
        if hasattr(self, "file") and self.file != None:
            self.file.close()
            self.file = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()



class BamAlignment(object):
    def __init__(self, bamReader, pysamAlignedRead):
        # Since BAM is a stream-oriented format, we load all the stuff
        # we need because we cannot (efficiently) seek back to the
        # record.
        dual             = pysamAlignedRead
        self.bam         = bamReader
        self.tId         = dual.tid
        self.tStart      = dual.pos
        self.tEnd        = dual.aend
        self.qStart      = dual.opt("XS")-1  # TODO: fix XS convention to be 0-based
        self.qEnd        = dual.opt("XE")
        self.qLen        = dual.qlen
        self.isReverse   = dual.is_reverse
        self.cigarString = dual.cigarstring  # Genome oriented
        self.query       = dual.query        # Genome oriented, excludes soft-clipped bases
        self.qName       = dual.qname
        self.MapQV       = dual.mapq

    @property
    def zmw(self):
        return self.zmwRead.zmw

    @property
    def zmwRead(self):
        if not self.cmpH5.moviesAttached:
            raise ValueError("Movies not attached!")
        return self.cmpH5.basH5Collection[self.readName]

    @property
    def rowNumber(self):
        # BAM/SAM doesn't really define this, which is a shame.  There
        # is a notion of "offset" in BAM, but it doesn't seem to be
        # exposed to pysam.
        return self.readName

    def clippedTo(self, refStart, refEnd):
        return ClippedBamAlignment(self, refStart, refEnd)

    @property
    def alignmentGroup(self):
        raise UnavailableFeature("BAM has no HDF5 groups")

    @property
    def referenceInfo(self):
        return self.bam.referenceInfo(self.RefGroupID)

    @property
    def referenceName(self):
        return self.referenceInfo.FullName

    @property
    def readName(self):
        # TODO: qname needs to be updated to have the s_e
        return "%s/%d_%d" % (self.qName, self.readStart, self.readEnd)

    @property
    def movieInfo(self):
        return self.cmpH5.movieInfo(self.MovieID)

    @property
    def isForwardStrand(self):
        return not self.isReverse

    @property
    def isReverseStrand(self):
        return self.isReverse

    @property
    def referenceId(self):
        return self.tId

    @property
    def referenceStart(self):
        return self.tStart

    @property
    def referenceEnd(self):
        return self.tEnd

    @property
    def readStart(self):
        return self.qStart

    @property
    def readEnd(self):
        return self.qEnd

    @property
    def referenceSpan(self):
        return self.tEnd - self.tStart

    @property
    def readLength(self):
        return self.qLen

    @property
    def alignedLength(self):
        raise Unimplemented()

    def spansReferencePosition(self, pos):
        return self.tStart <= pos < self.tEnd

    def spansReferenceRange(self, start, end):
        assert start <= end
        return (self.tStart <= start <= end <= self.tEnd)

    def overlapsReferenceRange(self, start, end):
        assert start <= end
        return (self.tStart < end) and (self.tEnd > start)

    def containedInReferenceRange(self, start, end):
        assert start <= end
        return (start <= self.tStart <= self.tEnd <= end)

    @property
    def accuracy(self):
        raise Unimplemented()

    @property
    def numPasses(self):
        raise Unimplemented()

    @property
    def zScore(self):
        raise UnavailableFeature("No ZScore in BAM")

    @property
    def barcode(self):
        raise Unimplemented()

    @property
    def barcodeName(self):
        raise Unimplemented()

    def alignmentArray(self, orientation="native"):
        raise UnavailableFeature()

    def transcript(self, orientation="native", style="gusfield"):
        raise Unimplemented()

    def cigar(self, orientation="native"):
        if orientation=="genomic" or self.isForwardStrand:
            return self.cigarString
        else:
            return reverseCigarString(self.cigarString)

    def read(self, aligned=True, orientation="native"):
        # SAM/BAM stores the query unaligned in genomic order
        shouldRC = self.isReverseStrand and orientation=="native"
        rawQuery = self.query
        if aligned:
            return readFromBamSeq(rawQuery, shouldRC, self.cigar(orientation="genomic"))
        else:
            return readFromBamSeq(rawQuery, shouldRC)

    def reference(self, aligned=True, orientation="native"):
        # BAM does not contain the actual reference, this returns "N"
        # with gaps inserted as appropriate
        seq = "N" * self.referenceSpan
        cigar = self.cigar(orientation="genomic") if aligned else None
        shouldRC = (orientation=="native") and (self.isReverseStrand)
        return refFromBamSeq(seq, shouldRC, cigar)

    def referencePositions(self, orientation="native"):
        # PYSAM BUG: as far as I can tell, "positions" is busted in pysam
        # We will have to manage this on our own.
        #referenceNonGapMask = (self.alignmentArray(orientation) & 0b1111) != GAP
        cigarUnpacked = np.fromstring(rlDecodeCigarString(self.cigar(orientation)),
                                      dtype=np.uint8)
        M = ord("M")
        D = ord("D")
        referenceNonGapMask = (cigarUnpacked == M) | (cigarUnpacked == D)
        if self.isReverseStrand and orientation == "native":
            return self.tEnd - 1 - np.hstack([0, np.cumsum(referenceNonGapMask[:-1])])
        else:
            return self.tStart + np.hstack([0, np.cumsum(referenceNonGapMask[:-1])])

    def readPositions(self, orientation="native"):
        raise Unimplemented()

    def pulseFeature(self, featureName, aligned=True, orientation="native"):
        raise Unimplemented()

    # IPD            = _makePulseFeatureAccessor("IPD")
    # PulseWidth     = _makePulseFeatureAccessor("PulseWidth")
    # QualityValue   = _makePulseFeatureAccessor("QualityValue")
    # InsertionQV    = _makePulseFeatureAccessor("InsertionQV")
    # DeletionQV     = _makePulseFeatureAccessor("DeletionQV")
    # DeletionTag    = _makePulseFeatureAccessor("DeletionTag")
    # MergeQV        = _makePulseFeatureAccessor("MergeQV")
    # SubstitutionQV = _makePulseFeatureAccessor("SubstitutionQV")

    # def __getattr__(self, key):
    #     return self.cmpH5.alignmentIndex[self.rowNumber][key]

    def __repr__(self):
        return "BAM alignment: %s  %3d  %9d  %9d" \
            % (("+" if self.isForwardStrand else "-"),
               self.referenceId, self.tStart, self.tEnd)

    def __str__(self): return repr(self)

    def __cmp__(self, other):
        return cmp((self.referenceId, self.tStart, self.tEnd),
                   (other.referenceId, other.tStart, other.tEnd))

    def __dir__(self):
        # Special magic improving IPython completion
        return ALIGNMENT_INDEX_COLUMNS


class ClippedBamAlignment(BamAlignment):

    def __init__(self, aln, winStart, winEnd):
        self.tId         = aln.tId
        self.MapQV       = aln.MapQV
        self.isReverse   = aln.isReverse
        self.qName       = aln.qName  # FIXME

        (self.query, self.qStart, self.qEnd,
         self.tStart, self.tEnd, self.cigarString) = computeClipping(aln, winStart, winEnd)

        self.qLen        = self.qEnd - self.qStart


# ------------------------------------------------------------
#  Helper functions.  All this stuff is really slow and should
#  almost certainly be moved into C unless we can find a very
#  clever way to do it with numpy

COMPLEMENT_MAP = { "A" : "T",
                   "T" : "A",
                   "C" : "G",
                   "G" : "C",
                   "N" : "N",
                   "-" : "-" }

def reverseComplement(seq):
    return "".join([ COMPLEMENT_MAP[b] for b in seq[::-1]])

def alignedRead(seq, cigar):
    alnChunks = []
    offset = 0
    for op in rlDecodeCigarString(cigar):
        if (op == "M") or (op == "I"):
            alnChunks.append(seq[offset:offset+1])
            offset += 1
        elif (op == "D"):
            alnChunks.append("-")
    return "".join(alnChunks)

def alignedRef(seq, cigar):
    alnChunks = []
    offset = 0
    for op in rlDecodeCigarString(cigar):
        if (op == "M") or (op == "D"):
            alnChunks.append(seq[offset:offset+1])
            offset += 1
        elif (op == "I"):
            alnChunks.append("-")
    return "".join(alnChunks)

def readFromBamSeq(seq, rc=False, cigar=None):
    if cigar is not None:
        seq = alignedRead(seq, cigar)
    if rc:
        seq = reverseComplement(seq)
    return seq

def refFromBamSeq(seq, rc=False, cigar=None):
    if cigar is not None:
        seq = alignedRef(seq, cigar)
    if rc:
        seq = reverseComplement(seq)
    return seq

DIGITS = map(str, range(10))

def rlDecodeCigarString(cigarStr):
    result = []
    runLength = 0
    for code in cigarStr:
        if code in DIGITS:
            n = int(code)
            if runLength == 0:
                runLength = n
            else:
                runLength = 10*runLength + n
        else:
            op = code
            result.append(op * runLength)
            runLength = 0
    return "".join(result)

def rlEncodeCigarString(cigarStr):
    return "".join("%d%s" % (len(list(ops)), op)
                   for (op, ops) in groupby(cigarStr))

def reverseCigarString(cigarStr):
    return rlEncodeCigarString(rlDecodeCigarString(cigarStr)[::-1])

def computeClipping(aln, winStart, winEnd):
    """
    For given reference clipping, returns
    (query, qStart, qEnd, tStart, tEnd, cigarString)
    """
    # Walk the cigar, update offset pointers;
    # Shadow parent class members
    refPos = aln.referenceStart
    readOffset = 0
    alnOffset = 0
    unpackedCigar = rlDecodeCigarString(aln.cigarString)

    # Extent pointers
    alnOffsetStart  = 0
    alnOffsetEnd    = len(unpackedCigar)
    readOffsetStart = 0
    readOffsetEnd   = aln.readLength

    for op in unpackedCigar:
        # Record start
        if refPos == winStart:
            alnOffsetStart = alnOffset
            readOffsetStart = readOffset
        # Update offsets
        alnOffset += 1
        assert op in "DIM"
        if op == "M":
            refPos += 1
            readOffset += 1
        elif op == "I":
            readOffset += 1
        elif op == "D":
            refPos += 1
        # Record end
        if refPos == winEnd:
            alnOffsetEnd = alnOffset
            readOffsetEnd = readOffset
            break

    offsetBegin = readOffsetStart
    offsetEnd   = readOffsetEnd
    qStart      = aln.qStart + readOffsetStart
    qEnd        = aln.qStart + readOffsetEnd
    tStart      = 5 # FIXME
    tEnd        = 6 # FIXME
    cigarString = rlEncodeCigarString(unpackedCigar[alnOffsetStart:alnOffsetEnd])
    query  = aln.query[offsetBegin:offsetEnd]

    return (query, qStart, qEnd, tStart, tEnd, cigarString)
