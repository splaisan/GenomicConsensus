
__all__ = ["BamReader", "BamAlignment"]

from pysam import Samfile
from pbcore.io import BasH5Collection, FastaTable
from pbcore.chemistry import decodeTriple, ChemistryLookupError


import numpy as np
from functools import wraps
from itertools import groupby
from os.path import abspath, expanduser
from bisect import bisect_right, bisect_left

class UnavailableFeature(Exception): pass
class Unimplemented(Exception):      pass
class ReferenceMismatch(Exception):  pass

PULSE_FEATURE_TAGS = { "InsertionQV"    : ("iq", "qv",   np.uint8),
                       "DeletionQV"     : ("dq", "qv",   np.uint8),
                       "DeletionTag"    : ("dt", "base", np.int8 ),
                       "SubstitutionQV" : ("sq", "qv",   np.uint8),
                       "MergeQV"        : ("mq", "qv",   np.uint8) }

class BamReader(object):

    def _loadReferenceInfo(self):
        refRecords = self.peer.header["SQ"]
        refNames   = [r["SN"] for r in refRecords]
        refLengths = [r["LN"] for r in refRecords]
        refMD5s    = [r["M5"] for r in refRecords]
        refIds = map(self.peer.gettid, refNames)
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


    def _loadReadGroupInfo(self):
        rgs = self.peer.header["RG"]
        readGroupTable_ = []
        pulseFeaturesInAll_ = frozenset(PULSE_FEATURE_TAGS.keys())
        for rg in rgs:
            rgID = rg["ID"]
            rgName = rg["PU"]
            ds = dict([pair.split("=") for pair in rg["DS"].split(";") if pair != ""])
            triple = ds["BINDINGKIT"], ds["SEQUENCINGKIT"], ds["SOFTWAREVERSION"]
            rgChem = decodeTriple(*triple)
            rgReadType = ds["READTYPE"]
            readGroupTable_.append((rgID, rgName, rgReadType, rgChem))
            pulseFeaturesInAll_ = pulseFeaturesInAll_.intersection(ds.keys())

        self._readGroupTable = np.rec.fromrecords(
            readGroupTable_,
            dtype=[("ID"                 , "O"),
                   ("MovieName"          , "O"),
                   ("ReadType"           , "O"),
                   ("SequencingChemistry", "O")])
        self._readGroupDict = { rg.ID : rg
                                for rg in self._readGroupTable }

        self._pulseFeaturesAvailable = pulseFeaturesInAll_


    def _loadProgramInfo(self):
        # TODO: guarantee that these fields are nonoptional in our bams --- check with Marcus
        # TODO: are we interesting in the PP info?
        self._programTable = np.rec.fromrecords(
            [ (pg["ID"], pg.get("VN", None), pg.get("CL", None))
              for pg in self.peer.header["PG"] ],
            dtype=[("ID"     ,     "O"),
                   ("Version",     "O"),
                   ("CommandLine", "O")])

    def _loadReferenceFasta(self, referenceFastaFname):
        ft = FastaTable(referenceFastaFname)
        # Verify that this FASTA is in agreement with the BAM's
        # reference table---BAM should be a subset.
        fastaIdsAndLens = set((c.id, c.length) for c in ft)
        bamIdsAndLens   = set((c.Name, c.Length) for c in self.referenceInfoTable)
        if not bamIdsAndLens.issubset(fastaIdsAndLens):
            raise ReferenceMismatch, "FASTA file must contain superset of reference contigs in BAM"
        self.referenceFasta = ft


    def __init__(self, fname, referenceFastaFname=None):
        self.filename = fname = abspath(expanduser(fname))
        self.peer = Samfile(fname, "rb")

        # Check for sortedness, index.
        # There doesn't seem to be a "public" way to do this right
        # now, but that's fine because we're going to have to rewrite
        # it all anyway once the pysam rewrite lands.
        if not self.peer._hasIndex:
            raise ValueError, "Specified bam file lacks a bam index---required for this API"

        self._loadReferenceInfo()
        self._loadReadGroupInfo()
        self._loadProgramInfo()

        if referenceFastaFname is not None:
            self._loadReferenceFasta(referenceFastaFname)


    @property
    def isIndexLoaded(self):
        return False

    @property
    def isReferenceLoaded(self):
        return self.referenceFasta is not None

    def attach(self, fofnFilename):
        self.basH5Collection = BasH5Collection(fofnFilename)

    @property
    def moviesAttached(self):
        return (self.basH5Collection is not None)

    @property
    def alignmentIndex(self):
        raise UnavailableFeature("BAM has no alignment index")

    #TODO: change concept to readGroupTable in cmp.h5
    @property
    def movieInfoTable(self):
        raise Unimplemented()

    # TODO: change to read group accessor, this is semantically wrong now
    def movieInfo(self, movieId):
        raise Unimplemented()

    @property
    def movieNames(self):
        return set([mi.MovieName for mi in self.readGroupTable])

    @property
    def readGroupTable(self):
        return self._readGroupTable

    def readGroup(self, readGroupId):
        return self._readGroupDict[readGroupId]

    @property
    def sequencingChemistry(self):
        """
        List of the sequencing chemistries by movie.  Order is
        unspecified.
        """
        return list(self.readGroupTable.SequencingChemistry)

    #TODO: elide "Info" innames?
    @property
    def referenceInfoTable(self):
        return self._referenceInfoTable

    #TODO: standard?  how about subread instead?  why capitalize ccs?
    # can we standardize this?  is cDNA an additional possibility
    @property
    def readType(self):
        """
        Either "standard", "CCS", "mixed", or "unknown", to represent the
        type of PacBio reads aligned in this BAM file.
        """
        readTypes = self.readGroupTable.ReadType
        if all(readTypes == "SUBREAD"):
            return "standard"
        elif all(readTypes == "CCS"):
            return "CCS"
        elif all((readTypes == "CCS") | (readTypes == "SUBREAD")):
            return "mixed"
        else:
            return "unknown"

    #TODO: Marcus needs to put something in the spec for this
    @property
    def version(self):
        raise Unimplemented()

    #TODO: Marcus needs to put something in the spec for this
    def versionAtLeast(self, minimalVersion):
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

    # TODO: make this private in cmp.h5 reader
    def alignmentGroup(self, alnGroupId):
        raise UnavailableFeature("BAM has no HDF5 groups")

    def referenceInfo(self, key):
        return self._referenceDict[key]

    # TODO: cmp.h5 readsInRange only accepts int key, not string.
    # that's just lame, fix it.
    def readsInRange(self, winId, winStart, winEnd, justIndices=False):
        # PYSAM BUG: fetch doesn't work if arg 1 is tid and not rname
        if not isinstance(winId, str):
            winId = self.peer.getrname(winId)
        if justIndices == True:
            raise UnavailableFeature("BAM is not random-access")
        else:
            return ( BamAlignment(self, it)
                     for it in self.peer.fetch(winId, winStart, winEnd, reopen=False) )

    def hasPulseFeature(self, featureName):
        return featureName in self._pulseFeaturesAvailable

    def pulseFeaturesAvailable(self):
        return self._pulseFeaturesAvailable

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
        self.peer.reset()
        for a in self.peer:
            yield BamAlignment(self, a)

    def __len__(self):
        return self.peer.mapped

    def close(self):
        if hasattr(self, "file") and self.file is not None:
            self.file.close()
            self.file = None

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.close()


def _makePulseFeatureAccessor(featureName):
    def f(self, aligned=True, orientation="native"):
        return self.pulseFeature(featureName, aligned, orientation)
    return f


def requiresReference(method):
    @wraps(method)
    def f(bamAln, *args, **kwargs):
        if not bamAln.bam.isReferenceLoaded:
            raise UnavailableFeature, "this feature requires loaded reference sequence"
        else:
            return method(bamAln, *args, **kwargs)
    return f

def requiresIndex(method):
    @wraps(method)
    def f(bamAln, *args, **kwargs):
        if not bamAln.bam.isIndexLoaded:
            raise UnavailableFeature, "this feature requires a PacBio BAM index"
        else:
            return method(bamAln, *args, **kwargs)
    return f

class BamAlignment(object):
    def __init__(self, bamReader, pysamAlignedRead):
        #TODO: make these __slot__
        self.peer        = pysamAlignedRead
        self.bam         = bamReader
        self.tStart      = self.peer.pos
        self.tEnd        = self.peer.aend
        # Our terminology doesn't agree with pysam's terminology for
        # "query", "read".  This makes this code confusing.
        if self.peer.is_reverse:
            clipLeft  = self.peer.rlen - self.peer.qend
            clipRight = self.peer.qstart
        else:
            clipLeft  = self.peer.qstart
            clipRight = self.peer.rlen - self.peer.qend
        self.rStart = self.qStart + clipLeft
        self.rEnd   = self.qEnd   - clipRight

    @property
    def qStart(self):
        return self.peer.opt("YS")

    @property
    def qEnd(self):
        return self.peer.opt("YE")

    @property
    def tId(self):
        return self.peer.tid

    @property
    def isReverseStrand(self):
        return self.peer.is_reverse

    @property
    def isForwardStrand(self):
        return not self.peer.is_reverse

    @property
    def HoleNumber(self):
        return self.peer.opt("ZM")

    @property
    def MapQV(self):
        return self.peer.mapq

    @property
    def zmw(self):
        return self.zmwRead.zmw

    @property
    def zmwRead(self):
        if not self.bam.moviesAttached:
            raise ValueError("Movies not attached!")
        return self.bam.basH5Collection[self.readName]

    # TODO: change name to "offset" to be generic
    @property
    def rowNumber(self):
        #raise Unimplemented()
        return "(unknown row)"

    def clippedTo(self, refStart, refEnd):
        """
        Return a new `BamAlignment` that refers to a subalignment of
        this alignment, as induced by clipping to reference
        coordinates `refStart` to `refEnd`.

        .. warning::
            This function takes time linear in the length of the alignment.
        """
        assert type(self) is BamAlignment
        if (refStart >= refEnd or
            refStart >= self.tEnd or
            refEnd   <= self.tStart):
            raise IndexError, "Clipping query does not overlap alignment"

        # The clipping region must intersect the alignment, though it
        # does not have to be contained wholly within it.
        refStart = max(self.referenceStart, refStart)
        refEnd   = min(self.referenceEnd,   refEnd)
        refPositions = self.referencePositions(orientation="genomic")
        readPositions = self.readPositions(orientation="genomic")
        uc = self.unrolledCigar(orientation="genomic")

        # Clipping positions within the alignment array
        clipStart = bisect_right(refPositions, refStart) - 1
        clipEnd   = bisect_left(refPositions, refEnd)

        # "The logic for setting rStart, rEnd is tragically
        # complicated, due to the end-exclusive coordinate system."
        tStart = refStart
        tEnd   = refEnd
        if self.isForwardStrand:
            rStart = readPositions[clipStart]
            rEnd   = readPositions[clipEnd - 1] + 1
            cUc = uc[clipStart:clipEnd]
        else:
            rStart = readPositions[clipEnd - 1]
            rEnd   = readPositions[clipStart] + 1
            cUc = uc[clipStart:clipEnd]
        return ClippedBamAlignment(self, tStart, tEnd, rStart, rEnd, cUc)

    #TODO: remove this
    @property
    def alignmentGroup(self):
        raise UnavailableFeature("BAM has no HDF5 groups")

    @property
    def referenceInfo(self):
        return self.bam.referenceInfo(self.referenceId)

    @property
    def referenceName(self):
        return self.referenceInfo.FullName

    @property
    def readName(self):
        if self.readGroup.ReadType == "CCS":
            return "%s/%d/%d_%d" % (self.readGroup.MovieName, self.HoleNumber, "ccs")
        else:
            return "%s/%d/%d_%d" % \
                (self.readGroup.MovieName, self.HoleNumber, self.readStart, self.readEnd)


    #TODO: get rid of this
    @property
    def movieInfo(self):
        raise Unimplemented()

    @property
    def readGroup(self):
        return self.bam.readGroup(self.peer.opt("RG"))

    @property
    def sequencingChemistry(self):
        return self.readGroup.SequencingChemistry

    @property
    def isForwardStrand(self):
        return not self.isReverseStrand

    @property
    def isReverseStrand(self):
        return self.peer.is_reverse

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
    def queryStart(self):
        return self.qStart

    @property
    def queryEnd(self):
        return self.qEnd

    #TODO: provide this in cmp.h5 but throw "unsupported"
    @property
    def queryName(self):
        return self.peer.qname

    @property
    def readStart(self):
        return self.rStart

    @property
    def readEnd(self):
        return self.rEnd

    @property
    def referenceSpan(self):
        return self.tEnd - self.tStart

    @property
    def readLength(self):
        return self.rEnd - self.rStart

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
    @requiresIndex
    def identity(self):
        if self.readLength == 0:
            return 0.
        else:
            return 1. - float(self.nMM + self.nIns + self.nDel)/self.readLength

    @property
    def numPasses(self):
        return self.peer.opt("NP")

    @property
    def zScore(self):
        raise UnavailableFeature("No ZScore in BAM")

    @property
    def barcode(self):
        raise Unimplemented()

    @property
    def barcodeName(self):
        raise Unimplemented()

    @requiresReference
    def transcript(self, orientation="native", style="gusfield"):
        """
        A text representation of the alignment moves (see Gusfield).
        This can be useful in pretty-printing an alignment.
        """
        uc = self.unrolledCigar(orientation)
        ref = np.fromstring(self.reference(aligned=True, orientation=orientation), dtype=np.int8)
        read = np.fromstring(self.read(aligned=True, orientation=orientation), dtype=np.int8)
        isMatch = (ref == read)

        # Disambiguate the "M" op
        cigarPlus = uc
        cigarPlus[(~isMatch) & (cigarPlus == BAM_CMATCH)] = BAM_CDIFF   # 'X'
        cigarPlus[( isMatch) & (cigarPlus == BAM_CMATCH)] = BAM_CEQUAL  # '='

        #                                    MIDNSHP=X
        _exoneratePlusTrans = np.fromstring("Z  ZZZZ|*", dtype=np.int8)
        _exonerateTrans     = np.fromstring("Z  ZZZZ| ", dtype=np.int8)
        _cigarTrans         = np.fromstring("ZIDZZZZMM", dtype=np.int8)
        _gusfieldTrans      = np.fromstring("ZIDZZZZMR", dtype=np.int8)

        if   style == "exonerate+": return _exoneratePlusTrans [cigarPlus].tostring()
        elif style == "exonerate":  return _exonerateTrans     [cigarPlus].tostring()
        elif style == "cigar":      return _cigarTrans         [cigarPlus].tostring()
        else:                       return _gusfieldTrans      [cigarPlus].tostring()


    @requiresReference
    def reference(self, aligned=True, orientation="native"):
        if not (orientation == "native" or orientation == "genomic"):
            raise ValueError, "Bad `orientation` value"
        tSeq = self.bam.referenceFasta[self.referenceName].sequence[self.tStart:self.tEnd]
        shouldRC = orientation == "native" and self.isReverseStrand
        tSeqOriented = reverseComplement(tSeq) if shouldRC else tSeq
        if aligned:
            x = np.fromstring(tSeqOriented, dtype=np.int8)
            y = self._gapifyRef(x, orientation)
            return y.tostring()
        else:
            return tSeqOriented

    def unrolledCigar(self, orientation="native"):
        """
        Run-length decode the CIGAR encoding, and orient.  Clipping ops are removed.
        """
        ucGenomic = unrollCigar(self.peer.cigar, exciseSoftClips=True)
        ucOriented = ucGenomic[::-1] if (orientation == "native" and self.isReverseStrand) else ucGenomic
        return ucOriented

    def referencePositions(self, aligned=True, orientation="native"):
        """
        Returns an array of reference positions.

        If aligned is True, the array has the same length as the
        alignment and referencePositions[i] = reference position of
        the i'th column in the oriented alignment.

        If aligned is False, the array has the same length as the read
        and referencePositions[i] = reference position of the i'th
        base in the oriented read.
        """
        assert (aligned in (True, False) and
                orientation in ("native", "genomic"))

        ucOriented = self.unrolledCigar(orientation)
        refNonGapMask = (ucOriented != BAM_CINS)

        if self.isReverseStrand and orientation == "native":
            pos = self.tEnd - 1 - np.hstack([0, np.cumsum(refNonGapMask[:-1])])
        else:
            pos = self.tStart + np.hstack([0, np.cumsum(refNonGapMask[:-1])])

        if aligned:
            return pos
        else:
            return pos[ucOriented != BAM_CDEL]

    def readPositions(self, aligned=True, orientation="native"):
        """
        Returns an array of read positions.

        If aligned is True, the array has the same length as the
        alignment and readPositions[i] = read position of the i'th
        column in the oriented alignment.

        If aligned is False, the array has the same length as the
        mapped reference segment and readPositions[i] = read position
        of the i'th base in the oriented reference segment.
        """
        assert (aligned in (True, False) and
                orientation in ("native", "genomic"))

        ucOriented = self.unrolledCigar(orientation)
        readNonGapMask = (ucOriented != BAM_CDEL)

        if self.isReverseStrand and orientation == "genomic":
            pos = self.rEnd - 1 - np.hstack([0, np.cumsum(readNonGapMask[:-1])])
        else:
            pos = self.rStart + np.hstack([0, np.cumsum(readNonGapMask[:-1])])

        if aligned:
            return pos
        else:
            return pos[ucOriented != BAM_CINS]


    def pulseFeature(self, featureName, aligned=True, orientation="native"):
        """
        Retrieve the pulse feature as indicated.
        - `aligned`    : whether gaps should be inserted to reflect the alignment
        - `orientation`: "native" or "genomic"
        """
        if not (orientation == "native" or orientation == "genomic"):
            raise ValueError, "Bad `orientation` value"
        if featureName == "read":
            kind_  = "base"
            dtype_ = np.int8
            data_  = self.peer.seq
        elif featureName == "QualityValue":
            kind_  = "raw"
            dtype_ = np.uint8
            data_  = self.peer.qual
        else:
            tag, kind_, dtype_ = PULSE_FEATURE_TAGS[featureName]
            data_ = self.peer.opt(tag)
        assert len(data_) == self.peer.rlen

        # In a SAM/BAM file, the read data is all reversed if the aln
        # is on the reverse strand.  Let's get it back in read
        # (native) orientation, and remove the other artifacts of BAM
        # encoding
        if self.isReverseStrand:
            if kind_ == "base": data = reverseComplement(data_)
            else:               data = data_[::-1]
        else:
            data = data_
        data = np.fromstring(data, dtype=dtype_)
        if kind_ == "qv": data -= 33
        del data_


        # [s, e) delimits the range, within the query, that is in the aligned read.
        # This will be determined by the soft clips actually in the file as well as those
        # imposed by the clipping API here.
        s = self.rStart - self.qStart
        e = self.rEnd   - self.qStart
        assert s >= 0 and e <= len(data)
        clipped = data[s:e]

        # How to present it to the user
        shouldReverse = self.isReverseStrand and orientation == "genomic"
        if kind_ == "base":
            ungapped = reverseComplementAscii(clipped) if shouldReverse else clipped
        else:
            ungapped = clipped[::-1] if shouldReverse else clipped

        if aligned == False:
            return ungapped
        else:
            return self._gapifyRead(ungapped, orientation)


    def _gapifyRead(self, data, orientation):
        return self._gapify(data, orientation, BAM_CDEL)

    def _gapifyRef(self, data, orientation):
        return self._gapify(data, orientation, BAM_CINS)

    def _gapify(self, data, orientation, gapOp):
        # Precondition: data must already be *in* the specified orientation
        if data.dtype == np.int8:
            gapCode = ord("-")
        else:
            gapCode = data.dtype.type(-1)
        uc = self.unrolledCigar(orientation=orientation)
        alnData = np.repeat(np.array(gapCode, dtype=data.dtype), len(uc))
        gapMask = (uc == gapOp)
        alnData[~gapMask] = data
        return alnData

    # TODO: We haven't yet decided where these guys are going to live.
    # IPD            = _makePulseFeatureAccessor("IPD")
    # PulseWidth     = _makePulseFeatureAccessor("PulseWidth")

    QualityValue   = _makePulseFeatureAccessor("QualityValue")
    InsertionQV    = _makePulseFeatureAccessor("InsertionQV")
    DeletionQV     = _makePulseFeatureAccessor("DeletionQV")
    DeletionTag    = _makePulseFeatureAccessor("DeletionTag")
    MergeQV        = _makePulseFeatureAccessor("MergeQV")
    SubstitutionQV = _makePulseFeatureAccessor("SubstitutionQV")

    def read(self, aligned=True, orientation="native"):
        feature = self.pulseFeature("read", aligned, orientation)
        return feature.tostring()

    def __repr__(self):
        return "BAM alignment: %s  %3d  %9d  %9d" \
            % (("+" if self.isForwardStrand else "-"),
               self.referenceId, self.tStart, self.tEnd)

    def __str__(self):
        if self.bam.isReferenceLoaded:
            COLUMNS = 80
            val = ""
            val += repr(self) + "\n\n"
            val += "Read:        " + self.readName           + "\n"
            val += "Reference:   " + self.referenceName      + "\n\n"
            val += "Read length: " + str(self.readLength)    + "\n"
            #val += "Identity:    " + "%0.3f" % self.identity + "\n"

            alignedRead = self.read()
            alignedRef = self.reference()
            transcript = self.transcript(style="exonerate+")
            refPos = self.referencePositions()
            refPosString = "".join([str(pos % 10) for pos in refPos])
            for i in xrange(0, len(alignedRef), COLUMNS):
                val += "\n"
                val += "  " + refPosString[i:i+COLUMNS] + "\n"
                val += "  " + alignedRef  [i:i+COLUMNS] + "\n"
                val += "  " + transcript  [i:i+COLUMNS] + "\n"
                val += "  " + alignedRead [i:i+COLUMNS] + "\n"
                val += "\n"
            return val
        else:
            return repr(self)

    def __cmp__(self, other):
        return cmp((self.referenceId, self.tStart, self.tEnd),
                   (other.referenceId, other.tStart, other.tEnd))


class ClippedBamAlignment(BamAlignment):
    def __init__(self, aln, tStart, tEnd, rStart, rEnd, unrolledCigar):

        # Self-consistency checks
        assert tStart <= tEnd
        assert rStart <= rEnd
        assert sum(unrolledCigar != BAM_CDEL) == (rEnd - rStart)

        self.peer   = aln.peer
        self.bam    = aln.bam
        self.tStart = tStart
        self.tEnd   = tEnd
        self.rStart = rStart
        self.rEnd   = rEnd
        self._unrolledCigar = unrolledCigar  # genomic orientation



    def unrolledCigar(self, orientation="native"):
        if orientation=="native" and self.isReverseStrand:
            return self._unrolledCigar[::-1]
        else:
            return self._unrolledCigar



# ------------------------------------------------------------
#  Helper functions.

COMPLEMENT_MAP = { "A" : "T",
                   "T" : "A",
                   "C" : "G",
                   "G" : "C",
                   "N" : "N",
                   "-" : "-" }

def complement(seq):
    return "".join([ COMPLEMENT_MAP[b] for b in seq ])

def reverseComplement(seq):
    return "".join([ COMPLEMENT_MAP[b] for b in seq[::-1]])

def complementAscii(a):
    return np.array([ord(COMPLEMENT_MAP[chr(b)]) for b in a], dtype=np.int8)

def reverseComplementAscii(a):
    return complementAscii(a)[::-1]


BAM_CMATCH     = 0
BAM_CINS       = 1
BAM_CDEL       = 2
BAM_CREF_SKIP  = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD       = 6
BAM_CEQUAL     = 7
BAM_CDIFF      = 8

def unrollCigar(cigar, exciseSoftClips=False):
    """
    Run-length decode the cigar (input is BAM packed CIGAR, not a cigar string)

    Removes hard clip ops from the output.  Remove all?
    """
    cigarArray = np.array(cigar, dtype=int)
    hasHardClipAtLeft = cigarArray[0,0] == BAM_CHARD_CLIP
    hasHardClipAtRight = cigarArray[-1,0] == BAM_CHARD_CLIP
    ncigar = len(cigarArray)
    x = np.s_[int(hasHardClipAtLeft) : ncigar - int(hasHardClipAtRight)]
    ops = np.repeat(cigarArray[x,0], cigarArray[x,1])
    if exciseSoftClips:
        return ops[ops != BAM_CSOFT_CLIP]
    else:
        return ops


# The routines below are unoptimized but are only ever used by code
# requesting data in an "aligned" context.  Quiver, at least, never
# asks for data like this.

def printAln(aln, fastaTable=None):
    query = aln.peer.seq
    uc = unrollCigar(aln.peer.cigarstring)
    tuples = []
    qpos = 0
    for op in uc:
        ref = "N"
        if op == BAM_CMATCH:
            tuples.append( (ref, query[qpos]) )
            qpos += 1
        elif op == BAM_CINS:
            tuples.append( ("-", query[qpos]) )
            qpos += 1
        elif op == BAM_CDEL:
            tuples.append( (ref, "-") )
        elif op == BAM_CSOFT_CLIP:
            qpos += 1
        else:
            raise Exception, "Unexpected CIGAR code"
    aref, aread = zip(*tuples)
    return "".join(aref), "".join(aread)
