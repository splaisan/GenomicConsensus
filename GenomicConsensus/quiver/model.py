import numpy as np
from GenomicConsensus.quiver.utils import asFloatFeature
import ConsensusCore as cc

__all__ = [ "ParameterSet",
            "AllQVsModel",
            "NoMergeQVModel",
            "NoQVsModel"      ]


class ParameterSet(object):
    def __init__(self, model, qvModelParams):
        self.model = model
        self.qvModelParams = qvModelParams

    @staticmethod
    def fromString(s):
        if   s == "NoQVsModel.trainedParams1":  return NoQVsModel.trainedParams1()
        elif s == "AllQVsModel.trainedParams1": return AllQVsModel.trainedParams1()
        elif s == "AllQVsModel.trainedParams2": return AllQVsModel.trainedParams2()
        elif s == "NoMergeQVModel.trainedParams1" : return NoMergeQVModel.trainedParams1()
        else: raise Exception, "Unrecognized parameter set"

    @staticmethod
    def bestAvailable(cmpH5):
        if AllQVsModel.isCompatibleWithCmpH5(cmpH5):
            params = AllQVsModel.trainedParams1()
        elif NoMergeQVModel.isCompatibleWithCmpH5(cmpH5):
            params = NoMergeQVModel.trainedParams1()
        else:
            params = NoQVsModel.trainedParams1()
        return params

class Model(object):
    @classmethod
    def paramsFromArray(cls, arr):
        assert len(arr) == cls.numFreeParams
        arr_ = np.zeros(shape=(14,))
        arr_[cls.freeParamIdx] = arr
        res_ = np.where(cls.fixedParamMask, cls.fullStart, arr_).astype(np.float32)
        return ParameterSet(cls, cc.QvModelParams(*res_.tolist()))

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
    def trainedParams1(cls):
        return cls.paramsFromArray(
            [ 10.51021999, -43.8755488,   -0.65519504,
              -24.11037889,  -1.07307557, -40.00499772,
              2.40005902,  -1.03174328,   -6.34582353,
              -1.76146179,  -40.9595257,   -4.854102 ])

    """
    Parameters from training against ref000001:10000-40000 @ 11x in
    job 038537, using the sigmoid objective function.
    """
    @classmethod
    def trainedParams2(cls):
        return cls.paramsFromArray(
            [2.6291355 , -27.33168616,  -0.39203815,
             -15.74840896, -0.55376003, -25.42894743,
             -10.50970384,  -0.67740848, -8.86229982,
             -1.66164741, -34.59014316,  -1.78042695])


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
    def trainedParams1(cls):
        return cls.paramsFromArray(
            [ -3.90937113e+01,  -2.52056402e+01,  -1.38876854e+00,
              -3.07886177e+01,  -1.51443235e-02,  -1.01846311e+00,
              -8.63561305e+00,  -1.86460591e+00,  -4.13471617e+01])


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
    def trainedParams1(cls):
        return cls.paramsFromArray(
            [-48.69212896,  -14.85421593,
              -10.00835906, -10.01483049,
              -14.85421593])
