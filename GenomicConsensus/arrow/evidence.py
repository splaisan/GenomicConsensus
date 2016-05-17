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

# Authors: David Alexander, Lance Hepler

__all__ = [ "dumpEvidence",
            "ArrowEvidence" ]

import h5py, logging, os.path, numpy as np
from collections import namedtuple
from itertools import groupby
from bisect import bisect_left, bisect_right
from pbcore.io import FastaReader, FastaWriter

import ConsensusCore2 as cc2

from .. import reference
from .utils import allSingleBaseMutations

_typeMap = { cc2.MutationType_INSERTION    : "Ins",
             cc2.MutationType_DELETION     : "Del",
             cc2.MutationType_SUBSTITUTION : "Sub" }

def _shortMutationDescription(mut, tpl):
    """
    More compact and uniform mutation description strings
    Examples:

    201 Ins . > G
    201 Sub C > T
    201 Del C > .
    """
    _type = _typeMap[mut.Type]
    _pos = mut.Start()
    _oldBase = "." if mut.Type == cc2.MutationType_INSERTION \
               else tpl[_pos]
    _newBase = "." if mut.Type == cc2.MutationType_DELETION \
               else mut.Base
    return "%d %s %s > %s" % (_pos, _type, _oldBase, _newBase)


def scoreMatrix(ai):
    """
    Returns (rowNames, columnNames, S)

    where:
      - S is a matrix where S_{ij} represents the score delta
        of mutation j against read i
      - rowNames[i] is an identifier name for the the read i---presently
        we use the the row number within the cmp.h5, encoded as a string
      - columnNames[j] is an identifier for mutation j, encoding the
        position, type, and base change
    """
    css = str(ai)
    allMutations = sorted(allSingleBaseMutations(css))
    readNames = list(ai.ReadNames())
    numReads = len(readNames)
    shape = (numReads, len(allMutations))
    scoreMatrix = np.zeros(shape)
    for j, mut in enumerate(allMutations):
        mutScores = ai.LLs(mut)
        scoreMatrix[:, j] = mutScores
    baselineScores =  np.array(ai.LLs())
    columnNames = [ _shortMutationDescription(mut, css)
                    for mut in allMutations ]
    rowNames = readNames
    return (rowNames, columnNames, baselineScores, scoreMatrix)


def dumpFullEvidence(evidenceDumpBaseDirectory,
                     refWindow, refSequence, alns,
                     arrowConsensus):
    # Format of evidence dump:
    # evidence_dump/
    #   ref000001/
    #     0-1005/
    #       reference.fa
    #       reads.fa
    #       consensus.fa
    #       arrow-scores.h5
    #     995-2005/
    #       ...
    join = os.path.join
    refId, refStart, refEnd = refWindow
    refName = reference.idToName(refId)
    windowDirectory = join(evidenceDumpBaseDirectory,
                           refName,
                           "%d-%d" % (refStart, refEnd))
    logging.info("Dumping evidence to %s" % (windowDirectory,))

    if os.path.exists(windowDirectory):
        raise Exception, "Evidence dump does not expect directory %s to exist." % windowDirectory
    os.makedirs(windowDirectory)
    refFasta       = FastaWriter(join(windowDirectory, "reference.fa"))
    readsFasta     = FastaWriter(join(windowDirectory, "reads.fa"))
    consensusFasta = FastaWriter(join(windowDirectory, "consensus.fa"))

    windowName = refName + (":%d-%d" % (refStart, refEnd))
    refFasta.writeRecord(windowName, refSequence)
    refFasta.close()

    consensusFasta.writeRecord(windowName + "|arrow", arrowConsensus.sequence)
    consensusFasta.close()

    rowNames, columnNames, baselineScores, scores = scoreMatrix(arrowConsensus.ai)
    arrowScoreFile = h5py.File(join(windowDirectory, "arrow-scores.h5"))
    arrowScoreFile.create_dataset("Scores", data=scores)
    vlen_str = h5py.special_dtype(vlen=str)
    arrowScoreFile.create_dataset("RowNames", data=rowNames, dtype=vlen_str)
    arrowScoreFile.create_dataset("ColumnNames", data=columnNames, dtype=vlen_str)
    arrowScoreFile.create_dataset("BaselineScores", data=baselineScores)
    arrowScoreFile.close()
    for aln in alns:
        readsFasta.writeRecord(str(aln.rowNumber),
                               aln.read(orientation="genomic", aligned=False))
    readsFasta.close()


def dumpFocusedEvidence(evidenceDumpBaseDirectory,
                        refWindow, refSequence, alns,
                        arrowConsensus, variants):

    join = os.path.join
    refId, refStart, refEnd = refWindow
    refName = reference.idToName(refId)
    windowDirectory = join(evidenceDumpBaseDirectory,
                           refName,
                           "%d-%d" % (refStart, refEnd))
    logging.info("Dumping evidence to %s" % (windowDirectory,))

    if os.path.exists(windowDirectory):
        raise Exception, "Evidence dump does not expect directory %s to exist." % windowDirectory
    os.makedirs(windowDirectory)

    # -- the work...

    ai = arrowConsensus.ai

    # Dump a more streamlined evidence with a focus on variants that
    # were called.

    if len(variants) > 1:
        logging.warning("Evidence dump: no support yet for multiple variants in a window.")
    variant = variants[0]

    def inverseMutation(refWindow, variant):
        # Identify the template mutation that "undoes" the variant
        if variant.variantType == "Insertion":
            return cc2.Mutation(cc2.MutationType_DELETION,
                                variant.refStart - refStart)
        else:
            raise NotImplementedError

    try:
        inv = inverseMutation(refWindow, variant)
    except NotImplementedError:
        logging.warn("Evidence dump: currently only support for insertion variants")
        return
    mutName = _shortMutationDescription(inv, arrowConsensus.sequence)

    scores = np.array(ai.LLs(inv))
    baselineScores = np.array(ai.LLs())
    readNames = ai.ReadNames()
    deltas = scores - baselineScores

    f = open(join(windowDirectory, "arrow-scores.csv"), "w")
    f.write("qName,aName,RowNumber,Read,Reference,ForwardOrientedRead," +
            "ForwardOrientedReference,ForwardAlignedRead,ForwardAlignedReference," +
            "SnrA,SnrC,SnrG,SnrT," +
            "BaselineLL," + "," + mutName + "\n")
    for (aln, readName, ll, llDelta) in zip(alns, readNames, baselineScores, deltas):
        assert readName == aln.readName  # FIXME: this will fail for inactive reads.
        identifiers = "%s,%s,%d" % (aln.qName, aln.readName, aln.rowNumber)
        bases = "%s,%s,%s,%s,%s,%s" % \
                (aln.read      (orientation="native",  aligned=False),
                 aln.reference (orientation="native",  aligned=False),
                 aln.read      (orientation="genomic", aligned=False),
                 aln.reference (orientation="genomic", aligned=False),
                 aln.read      (orientation="genomic", aligned=True),
                 aln.reference (orientation="genomic", aligned=True))
        snrMetrics = "%f,%f,%f,%f" % tuple(aln.hqRegionSnr)
        likelhoodMetrics = "%f,%f" % (ll, llDelta)

        f.write(",".join([identifiers,
                          bases,
                          snrMetrics,
                          likelhoodMetrics]) + "\n")

class ArrowEvidence(object):
    """
    An experimental reader class for arrow evidence dumps produced by
    arrow --dumpEvidence
    """

    Mutation = namedtuple("Mutation", ("Position", "Type", "FromBase", "ToBase"))

    @staticmethod
    def _parseMutName(mutName):
        fields = mutName.split(" ")
        pos = int(fields[0])
        type, fromBase, _, toBase = fields[1:]
        return ArrowEvidence.Mutation(pos, type, fromBase, toBase)

    def __init__(self, path, refStart, consensus, rowNames, colNames, baselineScores, scores):
        self.path           = path
        self.refStart       = refStart
        self.consensus      = consensus
        self.rowNames       = rowNames
        self.colNames       = colNames
        self.baselineScores = baselineScores
        self.scores         = scores
        self.muts           = map(ArrowEvidence._parseMutName, self.colNames)

    @property
    def positions(self):
        return  [ mut.Position for mut in self.muts ]

    @property
    def uniquePositions(self):
        return sorted(list(set(self.positions)))

    @property
    def delta(self):
        return self.scores - self.baselineScores[:, np.newaxis]

    @staticmethod
    def load(path):
        if path.endswith("/"): path = path[:-1]

        refWin_ = path.split("/")[-1].split("-")
        refStart = int(refWin_[0])

        with FastaReader(path + "/consensus.fa") as fr:
            consensus = next(iter(fr)).sequence

        with h5py.File(path + "/arrow-scores.h5", "r") as f:
            scores   = f["Scores"].value
            baselineScores = f["BaselineScores"].value
            colNames = f["ColumnNames"].value
            rowNames = f["RowNames"].value
            return ArrowEvidence(path, refStart, consensus,
                                 rowNames, colNames,
                                 baselineScores, scores)

    def forPosition(self, pos):
        posStart = bisect_left(self.positions, pos)
        posEnd   = bisect_right(self.positions, pos)
        return ArrowEvidence(self.path,
                             self.refStart,
                             self.consensus,
                             self.rowNames,
                             self.colNames[posStart:posEnd],
                             self.baselineScores,
                             self.scores[:, posStart:posEnd])


    def justSubstitutions(self):
        colMask = np.array(map(lambda s: ("Sub" in s), self.colNames))
        return ArrowEvidence(self.path,
                             self.refStart,
                             self.consensus,
                             self.rowNames,
                             self.colNames[colMask],
                             self.baselineScores,
                             self.scores[:, colMask])

    def rowNumbers(self):
        with FastaReader(self.path + "/reads.fa") as fr:
            return [ int(ctg.name) for ctg in fr ]
