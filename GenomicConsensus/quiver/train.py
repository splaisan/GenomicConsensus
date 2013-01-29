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

from multiprocessing import *
from ConsensusCore import (BandingOptions,
                           SparseSseQvMultiReadMutationScorer as MutationScorer)
from GenomicConsensus.quiver.utils import *
from GenomicConsensus.quiver.model import *

import argparse, math, numpy as np, pbls
from scipy.optimize import fmin, fmin_powell

TRAINING_JOB     = "038537"
TRAINING_WINDOWS = [ (1, winStart, winStart+1000)
                     for winStart in range(5000, 45000, 1000)]
BANDING   = BandingOptions(4, 200)
THRESHOLD = -500

def chunks(n, l):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def sigmoid(alternateScores):
    return sum(1./(1 + np.exp(-alternateScores)))

def logsigmoid(alternateScores):
    # Clip to prevent overflow in exp
    alternateScores_ = alternateScores = np.clip(alternateScores, -50, 50)
    return -sum(np.log(1 + np.exp(-alternateScores)))

def Worker_initialize():
    global cmpH5, refSeq
    cmpH5  = pbls.fetchCmpH5(TRAINING_JOB)
    refSeq = pbls.fetchReferenceSequence(cmpH5, "ref000001")

def Worker_objective(t):
    window, x = t
    loss = 0
    misclassifications = 0
    referenceBases = 0

    params = model.paramsFromArray(x, BANDING, THRESHOLD).quiverConfig

    if options.verbose: print "Window: " +  `window`
    _, winStart, winEnd = window
    trueTemplate = refSeq[winStart:winEnd]
    overlappingReads = cmpH5.readsInRange(*window)
    spanningReads = [aln for aln in overlappingReads
                     if aln.spansReferenceRange(*window[1:])]
    clippedReads = [aln.clippedTo(winStart, winEnd)
                    for aln in spanningReads]
    for usedReads in chunks(options.coverageDepth, clippedReads):
        # Skipping low-coverage regions
        if len(usedReads) < options.coverageDepth: continue

        mms = MutationScorer(params, trueTemplate)
        for aln in usedReads:
            qvsf = model.extractFeatures(aln)
            mms.AddRead(qvsf, int(aln.RCRefStrand))

        alternateScores = np.array([mms.Score(m)
                                    for m in uniqueSingleBaseMutations(trueTemplate)])

        # Apply a transformation and sum
        loss += transform(alternateScores)
        misclassifications += np.sum(alternateScores > 0)
        referenceBases += (winEnd - winStart)
    return (loss, misclassifications, referenceBases)

def initializeWorkers(numWorkers):
    return Pool(numWorkers, Worker_initialize)

fevals = 0

class memoize(object):
    """
    A very non-general memoize decorator class.  Works on a function
    of a single argument of type numpy.ndarray, which we have to
    pickle to a tuple (immutable, hashable) for storage in the cache.
    """
    def __init__(self,function):
        self.function=function
        self.cache = {}
    def __call__(self, arg):
        hashableArg = tuple(arg)
        if hashableArg in self.cache:
            print "Memoize cache hit!"
        else:
            self.cache[hashableArg] = self.function(arg)
        return self.cache[hashableArg]

# We memoize the objective function because scipy's optimizers seem to
# be pretty terrible at avoiding redundant function
# evaluations---which are very expensive in our case.

@memoize
def objective(x):
    # Process the data.
    results = pool.map(Worker_objective,
                       [(win, x) for win in TRAINING_WINDOWS])

    unzippedResults = zip(*results)
    loss, misclassifications, referenceBases = map(sum, unzippedResults)

    global fevals
    fevals += 1

    print "************"
    print "(Iteration #%d)"      % fevals
    print "x:                  " + str(x)
    print "Loss:               " + str(loss)
    print "Misclassifications: " + str(misclassifications)
    print "Quality:            " + \
        "Q" + str(int(-10*math.log10(float(misclassifications+1)/(referenceBases+1))))
    return loss

parser = argparse.ArgumentParser(description="Quiver training script")
parser.add_argument("--numWorkers", "-j",
                    action="store",
                    dest="numWorkers",
                    type=int,
                    default=8)
parser.add_argument("--coverageDepth", "-D",
                    action="store",
                    dest="coverageDepth",
                    type=int,
                    default=11)
parser.add_argument("--verbose", "-v",
                    action="store_true",
                    dest="verbose",
                    default=False)
parser.add_argument("--objective",
                    action="store",
                    dest="objectiveName",
                    default="logsigmoid",
                    choices=["logsigmoid",
                             "sigmoid"])
parser.add_argument("--model",
                    action="store",
                    dest="model",
                    default="AllQVsModel",
                    choices=["AllQVsModel",
                             "NoQVsModel",
                             "NoMergeQVModel"])

options = parser.parse_args()

transform = eval(options.objectiveName)
model = eval(options.model)

def main():
    np.seterr(over="ignore")

    global pool
    pool = initializeWorkers(options.numWorkers)
    x0 = model.start
    ans = fmin_powell(objective,
                      x0,
                      xtol=0.5,
                      full_output=True)
    print ans

if __name__=="__main__":
    main()
