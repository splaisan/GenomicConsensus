import cProfile, logging, os.path
from multiprocessing import Process
from threading import Thread
from pbcore.io import CmpH5Reader
from .options import options

class Worker(object):
    """
    Base class for compute worker that read reference coordinates
    from the task queue, perform variant calling, then push results
    back to another queue, to be written to a GFF file by a collector.

    All tasks that are O(genome length * coverage depth) should be
    distributed to compute workers, leaving the collector
    worker only O(genome length) work to do.
    """
    def __init__(self, workQueue, resultsQueue):
        self._workQueue = workQueue
        self._resultsQueue = resultsQueue

    def _run(self):
        self._inCmpH5 = CmpH5Reader(options.inputFilename)
        self.onStart()

        while True:
            datum = self._workQueue.get()
            if datum == None:
                # Sentinel indicating end of input.  Place a sentinel
                # on the results queue and end this worker process.
                self._resultsQueue.put(None)
                break
            else:
                coords, rowNumbers = datum
                logging.debug("%s received work unit, coords=%s, # reads=%d"
                              % (self.name, str(coords), len(rowNumbers)))

                alnHits = self._inCmpH5[rowNumbers]
                result = self.onChunk(coords, alnHits)
                self._resultsQueue.put(result)

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

    #==
    # Begin overridable interface
    #==

    def onStart(self):
        pass

    def onChunk(self, referenceWindow, alnHits):
        """
        This function is the heart of the matter.
        
        referenceWindow, alnHits -> result
        """
        pass

    def onFinish(self):
        pass

class WorkerProcess(Worker, Process):
    """Worker that executes as a process."""
    def __init__(self, *args):
        Process.__init__(self)
        super(WorkerProcess,self).__init__(*args)
        self.daemon = True

class WorkerThread(Worker, Thread):
    """Worker that executes as a thread (for debugging purposes only)."""
    def __init__(self, *args):
        Thread.__init__(self)
        super(WorkerThread,self).__init__(*args)
        self.daemon = True
        self.exitcode = 0
