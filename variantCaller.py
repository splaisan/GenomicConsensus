#!/usr/bin/env python

from __future__ import absolute_import
import argparse, atexit, cProfile, gc, glob, h5py, logging, multiprocessing
import os, pstats, random, shutil, sys, tempfile, time, threading, Queue
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

def die(msg):
    print msg
    sys.exit(-1)

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

        Raises an exception if file is unsorted or empty.
        """
        fname = options.inputFilename
        logging.info("Reading CmpH5 file %s" % (fname,))
        self._inCmpH5 = CmpH5Reader(fname)
        if not self._inCmpH5.isSorted:
            die("Input CmpH5 file must be sorted.")
        if self._inCmpH5.isEmpty:
            die("Input CmpH5 file must be nonempty.")

        # Compatible with selected algorithm?
        cmpH5isOK, msg = self._algorithm.compatibilityWithCmpH5(self._inCmpH5)
        if not cmpH5isOK:
            die("Failure: CmpH5 file is incompatible with algorithm \"%s\": %s" %
                (self._algorithm.name, msg))

        logging.info("Input CmpH5 data: numAlnHits=%d" % len(self._inCmpH5))

    def _loadReference(self):
        # The refGroupId<->refGroupFullName mapping is only delineated
        # in the cmph5 file.  So we have to open the cmph5, use it to
        # look up the mapping, then close it---we do the real open
        # after the fork.
        assert self._inCmpH5 == None
        cmph5 = CmpH5Reader(options.inputFilename)
        err = reference.loadFromFile(options.referenceFilename, cmph5)
        if err:
            die("Error loading reference")
        # Grok the referenceWindow spec, if any.
        options.referenceWindow = \
            reference.windowFromString(options.referenceWindowAsString)
        cmph5.close()

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

        self._loadReference()
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
                                filename=os.path.join(options.temporaryDirectory, "profile-main.out"))
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
