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

import csv, gzip, logging, os, numpy as np, sys, time
import traceback
from collections import namedtuple
from .. import reference
from ..variants import Insertion, Deletion, Substitution
from ..utils import fileFormat
from ..io  import VariantsGffWriter, FastqWriter, FastaWriter

# The output is written by "consumer" objects (opposite of generators)
# Each consumer consumes (refId, tbl) tuples, where the tbl is the numpy
# array containing the locus summary for each reference position.
# The only methods that should be needed here are:
#    - consumers
#    - broadcast

# ========================================
# Utility functions
#
def consumer(func):
    def start(*args,**kwargs):
        c = func(*args,**kwargs)
        c.next()
        return c
    start.__name__ = func.__name__
    start.__dict__ = func.__dict__
    start.__doc__  = func.__doc__
    return start

def broadcast(datum, consumers):
    for c in consumers:
        c.send(datum)

def printErrorMessage(where):
    head = "******** Error encountered in %s ********" % where
    print head
    traceback.print_exc()
    print "*" * len(head)
    print

# ========================================
# Consumers
#
@consumer
def csvConsumer(file, **kwargs):
    fieldNames = kwargs["fieldNames"]
    assert(fieldNames[0] == "referenceId")
    writer = csv.writer(file)
    writer.writerow(fieldNames)
    try:
        while True:
            (refId, tbl) = (yield)
            refHeader = reference.idToName(refId)
            coverage = tbl.view(np.recarray).coverage
            logging.info("Writing CSV output for %s." % refHeader)
            for record in tbl[coverage > 0, :]:
                writer.writerow((refHeader,) + record.item())
    except Exception:
        printErrorMessage("csvConsumer")
    finally:
        file.close()
        return

@consumer
def fastaConsumer(file, **kwargs):
    writer = FastaWriter(file)
    try:
        while True:
            (refId, tbl) = (yield)
            refHeader = reference.idToName(refId) + "|plurality"
            seqArray = tbl["consensus"]
            logging.info("Writing FASTA output for %s." % refHeader)
            writer.writeRecord(refHeader, seqArray)
    except Exception:
        printErrorMessage("fastaConsumer")
    finally:
        file.close()
        return

@consumer
def fastqConsumer(file, **kwargs):
    writer = FastqWriter(file)
    try:
        while True:
            (refId, tbl) = (yield)
            refHeader = reference.idToName(refId) + "|plurality"
            seqArray = tbl["consensus"]
            qvArray = tbl["consensusConfidence"]
            logging.info("Writing FASTQ output for %s." % refHeader)
            writer.writeRecord(refHeader, seqArray, qvArray)
    except Exception:
        printErrorMessage("fastqConsumer")
    finally:
        file.close()
        return

@consumer
def variantsGffConsumer(file, **kwargs):
    confidenceThreshold = kwargs["confidenceThreshold"]
    writer = VariantsGffWriter(file,
                               " ".join(sys.argv),
                               reference.byId.values())
    try:
        while True:
            (refId, tbl) = (yield)

            variants = []
            refHeader = reference.idToName(refId)
            refSeq = reference.byId[refId].sequence
            logging.info("Calculating variants from consensus for %s." % refHeader)

            loci = locusIterForGff(tbl, refSeq, confidenceThreshold)

            for locus in loci:
                if len(locus.consensus) == 0:
                    variants.append(Deletion(refId,
                                             locus.refPos,
                                             locus.refPos + 1,
                                             locus.refBase,
                                             "",
                                             coverage=locus.coverage,
                                             confidence=locus.confidence,
                                             frequency=locus.frequency))
                else:
                    if len(locus.consensus) > 1:
                        variants.append(Insertion(refId,
                                                  locus.refPos,
                                                  locus.refPos,
                                                  "",
                                                  locus.consensus[:-1],
                                                  coverage=locus.coverage,
                                                  confidence=locus.confidence,
                                                  frequency=locus.frequency))
                    if locus.consensus[-1] != locus.refBase:
                        variants.append(Substitution(refId,
                                                     locus.refPos,
                                                     locus.refPos+1,
                                                     locus.refBase,
                                                     locus.consensus[-1],
                                                     coverage=locus.coverage,
                                                     confidence=locus.confidence,
                                                     frequency=locus.frequency))

            logging.info("Writing variants for %s to GFF." % refHeader)
            for variant in variants:
                writer.writeRecord(variant.toGffRecord())

    except Exception:
        printErrorMessage("variantsGffConsumer")
    finally:
        file.close()
        return

Locus = namedtuple("Locus", ("refPos", "refBase",
                             "coverage", "consensus",
                             "confidence","frequency"))

def locusIterForGff(tbl, refSeq, qvThresh):
    """
    Obtain the approprate iterator for the given data structure.  The rare
    variant algorithm passes a dictionary.  All others are assumed to pass
    a numpy array (i.e., plurality, quiver).  It returns an iterator that
    returns a list of locus named tuples that can be used to generate the
    GFF file in a unified way.
    """
    # The rare variant algo passes a different data structure. Some
    # actual thought was put into trying to unify them, but it doesn't
    # make enough sense to spend too much time on it.
    if type(tbl) == dict:
        def _literDict(tbl, refSeq, qvThresh):
            for refPos, varl in tbl.iteritems():
                for v in varl:
                    yield Locus(refPos, refSeq[refPos], v.coverage, v.consensus,
                                  v.consensusConfidence, v.consensusFrequency)

        return _literDict(tbl, refSeq, qvThresh)

    else:
        def _literNumpy(tbl, refSeq, qvThresh):
            ra = tbl.view(np.recarray)
            # locates the variant positions in the reference
            varPos = np.flatnonzero((ra.consensus != refSeq) &
                                    (ra.consensusConfidence > qvThresh))
            for refPos in varPos:
                yield Locus(refPos, refSeq[refPos], ra.coverage[refPos],
                              ra.consensus[refPos], ra.consensusConfidence[refPos],
                              ra.consensusFrequency[refPos])

        return _literNumpy(tbl, refSeq, qvThresh)


def consumerForFilename(filename, **kwargs):
    if filename.endswith(".gz"):
        file = gzip.open(filename, "wb")
    else:
        file = open(filename, "w")

    format = fileFormat(filename)
    if format == "CSV":
        return csvConsumer(file, **kwargs)
    elif format == "FASTA":
        return fastaConsumer(file, **kwargs)
    elif format == "FASTQ":
        return fastqConsumer(file, **kwargs)
    elif format == "GFF":
        return variantsGffConsumer(file, **kwargs)
    else:
        raise ApplicationException, "Unsupported output format"


def consumers(filenames, fieldNames, confidenceThreshold=0):
    return [consumerForFilename(fn,
                                fieldNames=fieldNames,
                                confidenceThreshold=confidenceThreshold)
            for fn in filenames]


# ========================================
# Query for available output formats
#
def supportedOutputExtensions():
    uncompressedExtensions = [".csv",
                              ".fa", ".fasta",
                              ".fq", ".fastq",
                              ".gff"]
    compressedExtensions = [ext + ".gz" for ext in uncompressedExtensions]
    return uncompressedExtensions + compressedExtensions
