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
