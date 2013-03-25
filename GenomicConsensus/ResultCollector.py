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

# Author: David Alexander, Jim Drake

import cProfile, logging, os.path, sys
from multiprocessing import Process
from threading import Thread
from collections import OrderedDict, defaultdict
from .options import options
from GenomicConsensus import reference, consensus
from .io.VariantsGffWriter import VariantsGffWriter
from pbcore.io import FastaWriter, FastqWriter

class ResultCollector(object):
    """
    Gathers results and writes to a file.
    """
    def __init__(self, resultsQueue, algorithmConfig):
        self._resultsQueue = resultsQueue
        self._algorithmConfig = algorithmConfig

    def _run(self):
        self.onStart()

        sentinelsReceived = 0
        while sentinelsReceived < options.numWorkers:
            result = self._resultsQueue.get()
            if result is None:
                sentinelsReceived += 1
            else:
                self.onResult(result)

        self.onFinish()

    def run(self):
        if options.doProfiling:
            cProfile.runctx("self._run()",
                            globals=globals(),
                            locals=locals(),
                            filename=os.path.join(options.temporaryDirectory,
                                                  "profile-%s.out" % (self.name)))
        else:
            self._run()


    # ==================================
    # Overridable interface begins here.
    #

    def onStart(self):
        self.referenceBasesProcessedById = OrderedDict()
        for refId in reference.byId:
            self.referenceBasesProcessedById[refId] = 0
        self.variantsByRefId             = defaultdict(list)
        self.consensusChunksByRefId      = defaultdict(list)

        # open file writers
        self.fastaWriter = self.fastqWriter = self.gffWriter = None
        if options.fastaOutputFilename:
            self.fastaWriter = FastaWriter(options.fastaOutputFilename)
        if options.fastqOutputFilename:
            self.fastqWriter = FastqWriter(options.fastqOutputFilename)
        if options.gffOutputFilename:
            self.gffWriter = VariantsGffWriter(options.gffOutputFilename,
                                               " ".join(sys.argv),
                                               reference.byId.values())

    def onResult(self, result):
        window, cssAndVariants = result
        css, variants = cssAndVariants
        self._recordNewResults(window, css, variants)
        self._flushContigIfCompleted(window)

    def onFinish(self):
        logging.info("Analysis completed.")
        if self.fastaWriter: self.fastaWriter.close()
        if self.fastqWriter: self.fastqWriter.close()
        if self.gffWriter: self.gffWriter.close()
        logging.info("Output files completed.")

    def _recordNewResults(self, window, css, variants):
        refId, refStart, refEnd = window
        self.consensusChunksByRefId[refId].append(css)
        self.variantsByRefId[refId] += variants
        self.referenceBasesProcessedById[refId] += (refEnd - refStart)

    def _flushContigIfCompleted(self, window):
        refId, _, _ = window
        basesProcessed = self.referenceBasesProcessedById[refId]
        requiredBases = reference.numReferenceBases(refId, options.referenceWindow)
        if basesProcessed == requiredBases:
            # This contig is done, so we can dump to file and delete
            # the data structures.
            css = consensus.join(self.consensusChunksByRefId[refId])
            cssName = consensus.consensusContigName(reference.idToName(refId),
                                                    options.algorithm)
            if self.fastaWriter:
                self.fastaWriter.writeRecord(cssName,
                                             css.sequence)
            if self.fastqWriter:
                self.fastqWriter.writeRecord(cssName,
                                             css.sequence,
                                             css.confidence)
            del self.consensusChunksByRefId[refId]

            if self.gffWriter:
                self.gffWriter.writeVariants(sorted(self.variantsByRefId[refId]))
            del self.variantsByRefId[refId]

class ResultCollectorProcess(ResultCollector, Process):
    def __init__(self, *args):
        Process.__init__(self)
        self.daemon = True
        super(ResultCollectorProcess,self).__init__(*args)

class ResultCollectorThread(ResultCollector, Thread):
    def __init__(self, *args):
        Thread.__init__(self)
        self.daemon = True
        self.exitcode = 0
        super(ResultCollectorThread,self).__init__(*args)
