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
        return cmp( (self.refId, self.refStart, self.sortingPrecedence),
                    (other.refId, other.refStart, other.sortingPrecedence) )

    def __repr__(self):
        return "%d:%d-%d; %s -> %s" \
            % (self.refId, self.refStart, self.refEnd,
               self.refSequence or " ", self.readSequence or " ")

    def toTxtRecord(self):
        return repr(self)

    def toGffRecord(self):
        record = Gff3Record()
        record.seqid = reference.byId[self.refId].header
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
    sortingPrecedence = 1

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
    sortingPrecedence = 2

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
    sortingPrecedence = 3

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

