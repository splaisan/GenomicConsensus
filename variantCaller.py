#!/usr/bin/env python
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

__VERSION__ = "0.5.1"

#
# Verify that the correct version of the GenomicConsensus package is
# available before proceeding any further.
#
import sys
def die(msg):
    print msg
    sys.exit(-1)

try:
    from GenomicConsensus import __VERSION__ as packageVersion
    versionMismatch = (__VERSION__ != packageVersion)
except:
    versionMismatch = True
finally:
    if versionMismatch:
        die("\n"
            "Failure: The version of variantCaller.py does not match the version \n"
            "of the GenomicConsensus package.  The likely cause of this is       \n"
            "either:                                                             \n"
            "  1) You are running variantCaller.py without having activated the  \n"
            "     virtualenv containing the associated packages.                 \n"
            "     (Fix: activate the virtualenv)                                 \n"
            "  2) There is an older version of SMRTanalysis in your PYTHONPATH   \n"
            "     providing an older version of GenomicConsensus.                \n"
            "     (Fix: 'unset PYTHONPATH' in your shell)                        \n")

import argparse, atexit, cProfile, gc, glob, h5py, logging, multiprocessing
import os, pstats, random, shutil, tempfile, time, threading, Queue
from pbcore.io import CmpH5Reader

from GenomicConsensus import reference
from GenomicConsensus.options import (importAdditionalDefaultOptions,
                                      options,
                                      parseOptions,
                                      consensusCoreVersion)

from GenomicConsensus.quiver import quiver
from GenomicConsensus.plurality import plurality
from GenomicConsensus.rare import rare
from GenomicConsensus.utils import rowNumberIsInReadStratum

class ToolRunner(object):
    """
    The main driver class for the GenomicConsensus tool.
    """
    def __init__(self):
        self._inCmpH5 = None
        self._resultsQueue = None
        self._workQueue = None
        self._slaves = None
        self._algorithm = None
        self._aborting = False

    def _setupLogging(self):
        if options.verbosity >= 2:
            logLevel = logging.DEBUG
        elif options.verbosity == 1:
            logLevel = logging.INFO
        else:
            logLevel = logging.WARNING
        logFormat = '%(asctime)s [%(levelname)s] %(message)s'
        logging.basicConfig(level=logLevel, format=logFormat)

    def _makeTemporaryDirectory(self):
        """
        Make a temp dir where we can stash things if necessary.
        """
        options.temporaryDirectory = tempfile.mkdtemp(prefix="GenomicConsensus-", dir="/tmp")
        logging.info("Created temporary directory %s" % (options.temporaryDirectory,) )

    def _algorithmByName(self, name):
        if name=="plurality":
            algo = plurality
        elif name=="rare":
            algo = rare
        elif name=="quiver":
            algo = quiver
        else:
            die("Failure: unrecognized algorithm %s" % name)
        isOK, msg = algo.availability
        if not isOK:
            die("Failure: %s" % msg)
        return algo

    def _launchSlaves(self):
        """
        Launch a group of worker processes (self._slaves), the queue
        (self._workQueue) that will be used to send them chunks of
        work, and the queue that will be used to receive back the
        results (self._resultsQueue).

        Additionally, launch the result collector process.
        """
        availableCpus = multiprocessing.cpu_count()
        logging.info("Available CPUs: %d" % (availableCpus,))
        logging.info("Requested workers: %d" % (options.numWorkers,))
        logging.info("Parallel Mode: %s" % ("Threaded" if options.threaded else "Process",))
        if (options.numWorkers > availableCpus):
            logging.warn("More workers requested (%d) than CPUs available (%d);"
                         " may result in suboptimal performance."
                         % (options.numWorkers, availableCpus))
        self._initQueues()

        WorkerType, ResultCollectorType = self._algorithm.slaveFactories(options.threaded)
        self._slaves = []
        for i in xrange(options.numWorkers):
            p = WorkerType(self._workQueue, self._resultsQueue)
            self._slaves.append(p)
            p.start()
        logging.info("Launched compute slaves.")

        rcp = ResultCollectorType(self._resultsQueue)
        rcp.start()
        self._slaves.append(rcp)
        logging.info("Launched collector slave.")

    def _initQueues(self):
        if options.threaded:
            self._workQueue = Queue.Queue(options.queueSize)
            self._resultsQueue = Queue.Queue(options.queueSize)
        else:
            self._workQueue = multiprocessing.Queue(options.queueSize)
            self._resultsQueue = multiprocessing.Queue(options.queueSize)

    def _readCmpH5Input(self):
        """
        Read the CmpH5 input file into a CmpH5 object and
        store it as self._inCmpH5.
        """
        fname = options.inputFilename
        self._inCmpH5 = CmpH5Reader(fname)


    def _loadReference(self, cmpH5):
        err = reference.loadFromFile(options.referenceFilename, cmpH5)
        if err:
            die("Error loading reference")
        # Grok the referenceWindow spec, if any.
        options.referenceWindow = \
            reference.stringToWindow(options.referenceWindowAsString)

    def _checkFileCompatibility(self, cmpH5):
        if not cmpH5.isSorted:
            die("Input CmpH5 file must be sorted.")
        if cmpH5.isEmpty:
            die("Input CmpH5 file must be nonempty.")
        # Compatible with selected algorithm?
        cmpH5isOK, msg = self._algorithm.compatibilityWithCmpH5(cmpH5)
        if not cmpH5isOK:
            die("Failure: %s" %  msg)

    def _mainLoop(self):
        # Split up reference genome into chunks and farm out the
        # a chunk as a unit of work.
        logging.debug("Starting main loop.")
        ids = reference.enumerateIds(options.referenceWindow)
        for _id in ids:
            chunks = reference.enumerateChunks(_id,
                                               options.referenceChunkSize,
                                               options.referenceWindow,
                                               options.referenceChunkOverlap)
            for chunk in chunks:
                if self._aborting: return

                # Pickle row numbers, not alnHits
                rowNumbers = [rn for rn in self._inCmpH5.readsInRange(*chunk, justIndices=True)
                              if self._inCmpH5[rn].MapQV >= options.mapQvThreshold]

                if options.readStratum:
                    rowNumbers = [ rn for rn in rowNumbers
                                   if rowNumberIsInReadStratum(options.readStratum, rn) ]

                self._workQueue.put( (chunk, rowNumbers) )

        # Write sentinels ("end-of-work-stream")
        for i in xrange(options.numWorkers):
            self._workQueue.put(None)

    def _printProfiles(self):
        for profile in glob.glob(os.path.join(options.temporaryDirectory, "*")):
            pstats.Stats(profile).sort_stats("cumulative").print_stats(10)

    def _cleanup(self):
        if options.doProfiling:
            logging.info("Removing %s" % options.temporaryDirectory)
            shutil.rmtree(options.temporaryDirectory, ignore_errors=True)

    def _setupEvidenceDumpDirectory(self, directoryName):
        if os.path.exists(directoryName):
            shutil.rmtree(directoryName)
        os.makedirs(directoryName)

    @property
    def aborting(self):
        return self._aborting

    def abortWork(self):
        """
        Performs a shutdown of all the slave processes.  Called by the
        monitoring thread when a child process exits with a non-zero,
        or when a keyboard interrupt (Ctrl-C) is given. Not called
        during normal shutdown.
        """
        self._aborting = True
        self._resultsQueue.close()
        self._workQueue.close()

    @property
    def slaves(self):
        return self._slaves

    def main(self):
        # This looks scary but it's not.  Python uses reference
        # counting and has a secondary, optional garbage collector for
        # collecting garbage cycles.  Unfortunately when a cyclic GC
        # happens when a thread is calling cPickle.dumps, the
        # interpreter crashes sometimes.  See Bug 19704.  Since we
        # don't leak garbage cycles, disabling the cyclic GC is
        # essentially harmless.
        gc.disable()

        parseOptions()
        self._algorithm = self._algorithmByName(options.algorithm)
        importAdditionalDefaultOptions(self._algorithm.additionalDefaultOptions)

        self._setupLogging()
        random.seed(42)

        logging.info("h5py version: %s" % h5py.version.version)
        logging.info("hdf5 version: %s" % h5py.version.hdf5_version)
        logging.info("ConsensusCore version: %s" %
                     (consensusCoreVersion() or "ConsensusCore unavailable"))

        logging.info("Starting.")
        atexit.register(self._cleanup)
        if options.doProfiling:
            self._makeTemporaryDirectory()

        # We need to peek at the cmp.h5 file to build the The
        # refGroupId<->refGroupFullName mapping, and to determine
        # whether the selected algorithm parameters (Quiver) are
        # compatible with the data.  But we then have to close the
        # file, and let the "real" open happen after the fork.
        with CmpH5Reader(options.inputFilename) as peekCmpH5:
            logging.info("Peeking at CmpH5 file %s" % options.inputFilename)
            self._loadReference(peekCmpH5)
            self._checkFileCompatibility(peekCmpH5)
            logging.info("Input CmpH5 data: numAlnHits=%d" % len(peekCmpH5))

        if options.dumpEvidence:
            self._setupEvidenceDumpDirectory(options.evidenceDirectory)

        self._launchSlaves()
        self._readCmpH5Input()

        monitoringThread = threading.Thread(target=monitorSlaves, args=(self,))
        monitoringThread.start()

        try:
            if options.doProfiling:
                cProfile.runctx("self._mainLoop()",
                                globals=globals(),
                                locals=locals(),
                                filename=os.path.join(options.temporaryDirectory,
                                                      "profile-main.out"))
            else:
                self._mainLoop()
        except:
            self.abortWork()

        monitoringThread.join()

        if self._aborting:
            logging.error("Aborting")
            return -1
        else:
            logging.info("Finished.")

        if options.doProfiling:
            self._printProfiles()

        # close h5 file.
        del self._inCmpH5
        return 0


def monitorSlaves(driver):
    """
    Promptly aborts if a child is found to have exited with a nonzero
    exit code received; otherwise returns when all processes exit cleanly (0).

    This approach is portable--catching SIGCHLD doesn't work on
    Windows.
    """
    while not driver.aborting:
        all_exited = all(not p.is_alive() for p in driver.slaves)
        nonzero_exits = [p.exitcode for p in driver.slaves if p.exitcode]
        if nonzero_exits:
            exitcode = nonzero_exits[0]
            logging.error("Child process exited with exitcode=%d.  Aborting." % exitcode)
            driver.abortWork()
            return exitcode
        elif all_exited:
            return 0
        time.sleep(1)

if __name__ == "__main__":
    sys.exit(ToolRunner().main())
