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

from pbcore.io import Gff3Record, GffWriter
from . import reference
import sys, time

# Internally we use Python-style half-open intervals zero-based
# [start, end) to delineate reference ranges.  An insertion has
# start==end, a SNP has start+1==end, etc.
#
# GFF files assume 1-based indexing and open intervals [start, end).
# In a GFF both insertions and SNPs have start==end, which doesn't
# make too much sense to me, but so be it.

__all__ = [ "Deletion",
            "Insertion",
            "Substitution" ]

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
                 frequency=0,
                 zygosity=HAPLOID):
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
        self.coverage = coverage
        self.confidence = confidence
        self.frequency = frequency
        self.zygosity = zygosity

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
        gffStart, gffEnd = self.gffCoords()
        record = Gff3Record(reference.idToName(self.refId),
                            gffStart,
                            gffEnd,
                            self.gffType)
        if self.readSequence: record.variantSeq = self.readSequence
        if self.refSequence:  record.reference  = self.refSequence
        if self.zygosity != Variant.HAPLOID:
            record.zygosity = self.zygosity
        record.coverage = self.coverage
        record.confidence = self.confidence
        if self.frequency: record.frequency = self.frequency
        record.length = len(self)
        return record

class Insertion(Variant):
    sortingPrecedence = 1
    gffType = "insertion"

    def __repr__(self):
        return "Insertion   " + super(Insertion, self).__repr__()

    def __len__(self):
        return len(self.readSequence)

    def gffCoords(self):
        return self.refStart, self.refStart

class Deletion(Variant):
    sortingPrecedence = 2
    gffType = "deletion"

    def __repr__(self):
        return "Deletion    " + super(Deletion, self).__repr__()

    def __len__(self):
        return len(self.refSequence)

    def gffCoords(self):
        return self.refStart + 1, self.refEnd

class Substitution(Variant):
    sortingPrecedence = 3
    gffType = "substitution"

    def __repr__(self):
        return "Substitution" + super(Substitution, self).__repr__()

    def __len__(self):
        return len(self.refSequence)

    def gffCoords(self):
        return self.refStart + 1, self.refEnd
