from __future__ import absolute_import

import sys, time

from pbcore.io import Gff3Record, GffWriter
from . import reference

# Internally we use Python-style half-open intervals zero-based
# [start, end) to delineate reference ranges.  An insertion has
# start==end, a SNP has start+1==end, etc.
#
# GFF files assume 1-based indexing and open intervals [start, end).
# In a GFF both insertions and SNPs have start==end, which doesn't
# make too much sense to me, but so be it.

__all__ = [ "Deletion",
            "Insertion",
            "Substitution",
            "VariantList" ]

class Variant(object):

    HAPLOID = "haploid"
    HETEROZYGOUS = "heterozygous"
    HOMOZYGOUS = "homozygous"

    def __init__(self,
                 refId,
                 refStart, refEnd,
                 refSequence, readSequence,
                 coverage=0,
                 confidence=0,
                 zygosity=HAPLOID,
                 **kwargs):
        """
        `kwargs` encompasses the non-existential features of a
        variant, and are dumped as-is to GFF/CSV output file.  The
        variant caller must ensure that all variant calls have the
        same set of `kwargs` keys.
        """
        self.refId = refId
        self.refStart = refStart
        self.refEnd = refEnd
        self.refSequence = refSequence
        self.readSequence = readSequence
        self.zygosity = zygosity
        self.coverage = coverage
        self.confidence = confidence
        self.kwargs = kwargs

    def __cmp__(self, other):
        return cmp(self.refStart, other.refStart)

    def __repr__(self):
        return "%d:%d-%d; %s -> %s" \
            % (self.refId, self.refStart, self.refEnd,
               self.refSequence or " ", self.readSequence or " ")

    def toTxtRecord(self):
        return repr(self)

    def toGffRecord(self):
        record = Gff3Record()
        record.seqid = reference.byId[self.refId].name
        record.score = "."
        record.strand = "."

        # Empty attributes are not printed.
        if self.readSequence:
            record.put("variantSeq", self.readSequence)
        if self.refSequence:
            record.put("reference", self.refSequence)

        # If the zygosity attribute is omitted, HAPLOID is understood.
        if self.zygosity != Variant.HAPLOID:
            record.put("zygosity", self.zygosity)

        record.put("coverage", self.coverage)
        record.put("confidence", self.confidence)

        for key, value in self.kwargs.iteritems():
            record.put(key, value)
        return record


class Insertion(Variant):
    def __repr__(self):
        return "Insertion   " + super(Insertion, self).__repr__()

    def __len__(self):
        return len(self.readSequence)

    def toGffRecord(self):
        record = Variant.toGffRecord(self)
        record.type = "insertion"
        record.start, record.end = self.refStart, self.refStart
        record.put("length", len(self))
        return record

class Deletion(Variant):
    def __repr__(self):
        return "Deletion    " + super(Deletion, self).__repr__()

    def __len__(self):
        return len(self.refSequence)

    def toGffRecord(self):
        record = Variant.toGffRecord(self)
        record.type = "deletion"
        record.start, record.end = self.refStart + 1, self.refEnd
        record.put("length", len(self))
        return record

class Substitution(Variant):
    def __repr__(self):
        return "Substitution" + super(Substitution, self).__repr__()

    def __len__(self):
        return len(self.refSequence)

    def toGffRecord(self):
        record = Variant.toGffRecord(self)
        record.type = "substitution"
        record.start, record.end = self.refStart + 1, self.refEnd
        record.put("length", len(self))
        return record

class VariantList(object):
    """
    A collection class for variants that allows
    sorting and basic cleanups based on compacting
    adjacent variants when possible.
    """
    def __init__(self):
        self._lst = []

    def __iter__(self):
        return iter(self._lst)

    def __repr__(self):
        return "\n".join(repr(v) for v in self)

    def add(self, variant):
        self._lst.append(variant)

    def addDeletion(self, refId, refStart, refEnd, **kwargs):
        refSequence = reference.byId[refId].sequence[refStart:refEnd].tostring()
        self.add(Deletion(refId, refStart, refEnd, refSequence, "", **kwargs))

    def addSubstitution(self, refId, refStart, refEnd, readSequence, **kwargs):
        assert len(readSequence) == refEnd - refStart
        refSequence = reference.byId[refId].sequence[refStart:refEnd].tostring()
        self.add(Substitution(refId, refStart, refEnd, refSequence, readSequence, **kwargs))

    def addInsertion(self, refId, refStart, refEnd, readSequence, **kwargs):
        assert refStart == refEnd
        self.add(Insertion(refId, refStart, refEnd, "", readSequence, **kwargs))

    def sort(self):
        self._lst.sort()

