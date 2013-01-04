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

import numpy as np
from GenomicConsensus.quiver.utils import asFloatFeature
import ConsensusCore as cc

__all__ = [ "ParameterSet",
            "AllQVsModel",
            "NoMergeQVModel",
            "NoQVsModel"      ]


class ParameterSet(object):
    def __init__(self, model, quiverConfig):
        self.model = model
        self.quiverConfig = quiverConfig

    @staticmethod
    def fromString(s):
        if   s == "NoQVsModel.C2":      return NoQVsModel.C2()
        elif s == "AllQVsModel.C2":     return AllQVsModel.C2()
        elif s == "NoMergeQVModel.C2" : return NoMergeQVModel.C2()
        else: raise Exception, "Unrecognized parameter set"

    @staticmethod
    def bestAvailable(cmpH5):
        if AllQVsModel.isCompatibleWithCmpH5(cmpH5):
            params = AllQVsModel.C2()
        elif NoMergeQVModel.isCompatibleWithCmpH5(cmpH5):
            params = NoMergeQVModel.C2()
        else:
            params = NoQVsModel.C2()
        return params

class Model(object):
    @classmethod
    def paramsFromArray(cls, arr, bandingOptions, fastScoreThreshold):
        assert len(arr) == cls.numFreeParams
        arr_ = np.zeros(shape=(14,))
        arr_[cls.freeParamIdx] = arr
        res_ = np.where(cls.fixedParamMask, cls.fullStart, arr_).astype(np.float32)
        qvModelParams = cc.QvModelParams(*res_.tolist())
        return ParameterSet(cls, cc.QuiverConfig(qvModelParams,
                                                 cc.ALL_MOVES,
                                                 bandingOptions,
                                                 fastScoreThreshold))

    requiredFeatures = set([])

    @classmethod
    def isCompatibleWithCmpH5(cls, cmpH5):
        return all(cmpH5.hasPulseFeature(feature) for feature in cls.requiredFeatures)

    @classmethod
    def extractFeatures(cls, aln):
        """
        Extract the data in a cmp.h5 alignment record into a
        ConsensusCore-friendly `QvSequenceFeatures` object.  Will
        extract only the features relevant to this Model, zero-filling
        the other features arrays.

        Note that we have to use the AlnArray to see where the gaps
        are, at least for the moment (see bug 20752).
        """
        alnRead = np.fromstring(aln.read(), dtype=np.int8)
        gapMask = alnRead == ord("-")
        _args = [ alnRead[~gapMask].tostring() ]
        for feature in [ "InsertionQV",
                         "SubstitutionQV",
                         "DeletionQV",
                         "DeletionTag",
                         "MergeQV" ]:
            if feature in cls.requiredFeatures:
                _args.append(asFloatFeature(aln.pulseFeature(feature)[~gapMask]))
            else:
                _args.append(cc.FloatFeature(int(aln.readLength)))
        return cc.QvSequenceFeatures(*_args)

    @classmethod
    def extractMappedRead(cls, aln, windowStart):
        """
        Given a clipped alignment, convert its coordinates into template
        space (starts with 0), bundle it up with its features as a
        MappedRead.
        """
        assert aln.referenceSpan > 0
        return cc.MappedRead(cls.extractFeatures(aln),
                             int(aln.RCRefStrand),
                             int(aln.referenceStart) - windowStart,
                             int(aln.referenceEnd)   - windowStart)


class AllQVsModel(Model):

    requiredFeatures = set([ "InsertionQV",
                             "SubstitutionQV",
                             "DeletionQV",
                             "DeletionTag",
                             "MergeQV"       ])

    freeParamIdx   = range(12)  # Everything but the Burst stuff
    fixedParamIdx  = []
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(14) ]
    numFreeParams  = len(freeParamIdx)
    fullStart      = np.array([ -1.28069103,  -31.42568016,
                                 -0.44760609, -15.58668041,
                                 -0.4832615 , -25.06574059,
                                 -3.10107303,  -0.71846718,
                                 -18.,         -1.70749295,
                                 -44.3414917,   0.,
                                 0.,            0.         ],
                              dtype=np.float32)

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]

    """
    Parameters from training against ref000001:10000-40000 @ 11x in
    job 038537, using the logsigmoid objective function.
    """
    @classmethod
    def C2(cls):
        return cls.paramsFromArray(
            np.array([ 0.2627555 , -1.09688872, -0.01637988, -0.60275947, -0.02682689,
                       -1.00012494,  0.06000148, -0.02579358, -0.15864559, -0.04403654,
                       -1.02398814, -0.12135255]),
            bandingOptions=cc.BandingOptions(4, 5),
            fastScoreThreshold=-12.5)


class NoMergeQVModel(Model):
    """
    This model is intended for cmp.h5 files produced using the
    ResequencingQVs workflow using bas.h5 files that lack the MergeQV
    (i.e. Primary software pre-1.3.1).
    """
    requiredFeatures = set([ "InsertionQV",
                             "SubstitutionQV",
                             "DeletionQV",
                             "DeletionTag"])

    freeParamIdx = [ 1,  # Mismatch
                     3,  # Branch
                     4,  # BranchS
                     5,  # DeletionN
                     6,  # DeletionWithTag
                     7,  # DeletionWithTagS
                     8,  # Nce
                     9,  # NceS
                    10 ] # Merge

    fixedParamIdx = [ i for i in xrange(14) if i not in freeParamIdx ]
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(14) ]
    numFreeParams = len(freeParamIdx)
    fullStart = AllQVsModel.fullStart

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]

    """
    Parameters from training against ref000001:10000-40000 @ 11x in
    job 038537, using the logsigmoid objective function.
    """
    @classmethod
    def C2(cls):
        return cls.paramsFromArray(
            [ -3.90937113e+01,  -2.52056402e+01,  -1.38876854e+00,
              -3.07886177e+01,  -1.51443235e-02,  -1.01846311e+00,
              -8.63561305e+00,  -1.86460591e+00,  -4.13471617e+01],
            bandingOptions=cc.BandingOptions(4, 200),
            fastScoreThreshold=-500)


class NoQVsModel(Model):

    requiredFeatures = set([])

    freeParamIdx =   [ 1,   # Mismatch
                       3,   # Branch;
                       5,   # DeletionN;
                       8,   # Nce;
                       10 ] # Merge;

    fixedParamIdx = [ i for i in xrange(14) if i not in freeParamIdx ]
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(14) ]
    numFreeParams = len(freeParamIdx)

    fullStart = -10*np.array(~np.array(fixedParamMask), dtype=np.float32)

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]

    """
    Parameters from training against ref000001:10000-40000 @ 11x in
    job 038537, using the logsigmoid objective function.
    """
    @classmethod
    def C2(cls):
        return cls.paramsFromArray(
            [-48.69212896,  -14.85421593,
              -10.00835906, -10.01483049,
              -14.85421593],
            bandingOptions=cc.BandingOptions(4, 200),
            fastScoreThreshold=-500)
