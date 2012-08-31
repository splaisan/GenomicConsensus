import cProfile, logging, os.path
from multiprocessing import Process
from threading import Thread
from .options import options

class ResultCollector(object):
    """
    Gathers results and writes to a file.
    """
    def __init__(self, resultsQueue):
        self._resultsQueue = resultsQueue

    def _run(self):
        self.onStart()

        sentinelsReceived = 0
        while sentinelsReceived < options.numWorkers:
            result = self._resultsQueue.get()
            if result == None:
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
        pass

    def onResult(self, result):
        pass

    def onFinish(self):
        pass

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
