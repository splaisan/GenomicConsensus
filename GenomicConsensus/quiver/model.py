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

from GenomicConsensus.utils import die
from GenomicConsensus.quiver.utils import asFloatFeature, fst, snd
import ConsensusCore as cc

import numpy as np, ConfigParser, collections
from glob import glob
from os.path import join
from pkg_resources import resource_filename, Requirement


__all__ = [ "ParameterSet",
            "AllQVsModel",
            "NoMergeQVModel",
            "NoQVsModel",
            "AllQVsMergingByChannelModel",
            "NoQVsMergingByChannelModel",
            "findParametersFile",
            "loadParameterSets",
            "bestParameterSet",
            "majorityChemistry" ]


_basicParameterNames = [ "Match",
                         "Mismatch",
                         "MismatchS",
                         "Branch",
                         "BranchS",
                         "DeletionN",
                         "DeletionWithTag",
                         "DeletionWithTagS",
                         "Nce",
                         "NceS",
                         "Merge",
                         "MergeS" ]

_mergeByChannelParameterNames = [ "Match",
                                  "Mismatch",
                                  "MismatchS",
                                  "Branch",
                                  "BranchS",
                                  "DeletionN",
                                  "DeletionWithTag",
                                  "DeletionWithTagS",
                                  "Nce",
                                  "NceS",
                                  "Merge_A",
                                  "Merge_C",
                                  "Merge_G",
                                  "Merge_T",
                                  "MergeS_A",
                                  "MergeS_C",
                                  "MergeS_G",
                                  "MergeS_T" ]


class ParameterSet(object):
    def __init__(self, name, model, chemistry, quiverConfig):
        self.name         = name
        self.chemistry    = chemistry
        self.model        = model
        self.quiverConfig = quiverConfig

def _getResourcesDirectory():
    return resource_filename(Requirement.parse("GenomicConsensus"),
                             "GenomicConsensus/quiver/resources")

def majorityChemistry(cmpH5):
    """
    For the moment, we are doing Quiver analyses based on the majority
    chemistry represented in the cmp.h5 file.  Admittedly this could
    lead to suboptimal Quiver performance on mixed-chemistry cmp.h5
    files, but it is expedient.  Tie-breaking is done by alphabetical
    order of chemistry name.
    """
    chemistries = cmpH5.movieInfoTable.SequencingChemistry
    counts = collections.Counter(chemistries).most_common()
    sortedCounts = sorted(counts, key=lambda t: (t[1], t[0]), reverse=True)
    return sortedCounts[0][0]

def findParametersFile(filenameOrDirectory=None):
    if filenameOrDirectory is None:
        filenameOrDirectory = _getResourcesDirectory()

    # Given a full path to an .ini file, return the path
    if filenameOrDirectory.endswith(".ini"):
        return filenameOrDirectory

    # Given a path to a bundle (the directory with a date as its
    # name), return the path to the .ini file within
    foundInThisBundle = glob(join(filenameOrDirectory,
                                  "GenomicConsensus/QuiverParameters.ini"))
    if foundInThisBundle:
        return foundInThisBundle[0]

    # Given a directory containing bundles, return the path to the
    # .ini file within the lexically largest bundle subdirectory
    foundInBundlesBelow = glob(join(filenameOrDirectory,
                                    "*/GenomicConsensus/QuiverParameters.ini"))
    if foundInBundlesBelow:
        return sorted(foundInBundlesBelow)[-1]

    raise ValueError("Unable to find parameter set file (QuiverParameters.ini)")

def _buildParameterSet(parameterSetName, nameValuePairs):
    chem, modelName = parameterSetName.split(".")[:2]
    if    modelName == "AllQVsModel":    model = AllQVsModel
    elif  modelName == "NoMergeQVModel": model = NoMergeQVModel
    elif  modelName == "NoQVsModel":     model = NoQVsModel
    elif  modelName == "AllQVsMergingByChannelModel": model = AllQVsMergingByChannelModel
    elif  modelName == "NoQVsMergingByChannelModel":  model = NoQVsMergingByChannelModel
    else:
        logging.error("Found parameter set for unrecognized model: %s" % modelName)
        return None

    if map(fst, nameValuePairs) != model.parameterNames:
        die("Malformed parameter set file")

    qvModelParams = cc.QvModelParams(*[ float(snd(pair)) for pair in nameValuePairs ])
    quiverConfig = cc.QuiverConfig(qvModelParams,
                                   cc.ALL_MOVES,
                                   cc.BandingOptions(4, 5),
                                   -12.5)
    return ParameterSet(parameterSetName, model, chem, quiverConfig)

def loadParameterSets(iniFilename):
    # returns dict: name -> ParameterSet
    cp = ConfigParser.ConfigParser()
    cp.optionxform=str
    cp.read([iniFilename])
    sections = cp.sections()
    parameterSets = {}
    for sectionName in sections:
        parameterSet = _buildParameterSet(sectionName, cp.items(sectionName))
        if parameterSet:
            parameterSets[sectionName] = parameterSet
    return parameterSets

def bestParameterSet(parameterSets, chemistry, qvsAvailable):
    fallbackParameterSets = \
        [ paramSet for paramSet in parameterSets
          if paramSet.chemistry == "unknown"
          if paramSet.model.requiredFeatures.issubset(qvsAvailable) ]
    perChemistryParameterSets = \
        [ paramSet for paramSet in parameterSets
          if paramSet.chemistry == chemistry
          if paramSet.model.requiredFeatures.issubset(qvsAvailable) ]
    # Find the best one, under the assumption that a chemistry-trained
    # parameter set is always better than the "unknown" chemistry set.
    if perChemistryParameterSets:
        return max(perChemistryParameterSets, key=lambda ps: ps.model.rank)
    elif fallbackParameterSets:
        return max(fallbackParameterSets,     key=lambda ps: ps.model.rank)
    else:
        raise Exception("Quiver: No applicable parameter set found!")


class Model(object):
    @classmethod
    def paramsFromArray(cls, arr, bandingOptions, fastScoreThreshold):
        assert len(arr) == cls.numFreeParams
        arr_ = np.zeros(shape=(len(cls.parameterNames),))
        arr_[cls.freeParamIdx] = arr
        res_ = np.where(cls.fixedParamMask, cls.fullStart, arr_).astype(np.float32)
        qvModelParams = cc.QvModelParams(*res_.tolist())
        return ParameterSet("Training", "Training", cls,
                            cc.QuiverConfig(qvModelParams,
                                            cc.ALL_MOVES,
                                            bandingOptions,
                                            fastScoreThreshold))

    requiredFeatures = set([])
    parameterNames = []

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
    # Rank is used to determine whether one model is better than another,
    # all else being equal
    rank = 3

    requiredFeatures = set([ "InsertionQV",
                             "SubstitutionQV",
                             "DeletionQV",
                             "DeletionTag",
                             "MergeQV"       ])

    parameterNames = _basicParameterNames

    freeParamIdx   = range(12)
    fixedParamIdx  = []
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(12) ]
    numFreeParams  = len(freeParamIdx)

    #
    # This is the C2 parameter set, which will also be used as the
    # starting point for training.  These parameters are from training
    # against ref000001:10000-40000 @ 11x in job 038537, using the
    # logsigmoid objective function.
    #
    fullStart = np.array([ 0.2627555 , -1.09688872,
                           -0.01637988, -0.60275947,
                           -0.02682689, -1.00012494,
                           0.06000148, -0.02579358,
                           -0.15864559, -0.04403654,
                           -1.02398814, -0.12135255],
                         dtype=np.float)

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]

class NoMergeQVModel(Model):
    """
    This model is intended for cmp.h5 files produced using the
    ResequencingQVs workflow using bas.h5 files that lack the MergeQV
    (i.e. Primary software pre-1.3.1).
    """
    rank = 2

    requiredFeatures = set([ "InsertionQV",
                             "SubstitutionQV",
                             "DeletionQV",
                             "DeletionTag"])

    parameterNames = _basicParameterNames

    freeParamIdx = [ 0,  # Match
                     1,  # Mismatch
                     2,  # MismatchS
                     3,  # Branch
                     4,  # BranchS
                     5,  # DeletionN
                     6,  # DeletionWithTag
                     7,  # DeletionWithTagS
                     8,  # Nce
                     9,  # NceS
                    10 ] # Merge

    fixedParamIdx = [ i for i in xrange(12) if i not in freeParamIdx ]
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(12) ]
    numFreeParams = len(freeParamIdx)
    fullStart = AllQVsModel.fullStart

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]


class NoQVsModel(Model):
    rank = 1
    requiredFeatures = set([])
    parameterNames = _basicParameterNames

    freeParamIdx =   [ 1,   # Mismatch
                       3,   # Branch;
                       5,   # DeletionN;
                       8,   # Nce;
                       10 ] # Merge;

    fixedParamIdx = [ i for i in xrange(12) if i not in freeParamIdx ]
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(12) ]
    numFreeParams = len(freeParamIdx)

    fullStart = -10*np.array(~np.array(fixedParamMask), dtype=np.float32)

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]


class AllQVsMergingByChannelModel(Model):
    rank = -1
    requiredFeatures = set([ "InsertionQV",
                             "SubstitutionQV",
                             "DeletionQV",
                             "DeletionTag",
                             "MergeQV"       ])

    parameterNames = _mergeByChannelParameterNames


    freeParamIdx   = range(18)
    fixedParamIdx  = []
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(18) ]
    numFreeParams  = len(freeParamIdx)

    #
    # This is the C2 parameter set, which will also be used as the
    # starting point for training.  These parameters are from training
    # against ref000001:10000-40000 @ 11x in job 038537, using the
    # logsigmoid objective function.
    #
    fullStart = np.array([ 0.2627555 , -1.09688872,
                           -0.01637988, -0.60275947,
                           -0.02682689, -1.00012494,
                           0.06000148, -0.02579358,
                           -0.15864559, -0.04403654,
                           -1.02398814, -1.02398814,
                           -1.02398814, -1.02398814,
                           -0.12135255, -0.12135255,
                           -0.12135255, -0.12135255 ],
                         dtype=np.float)

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]

class NoQVsMergingByChannelModel(Model):
    rank = -1
    requiredFeatures = set([])
    parameterNames = _mergeByChannelParameterNames

    freeParamIdx =   [ 1,   # Mismatch
                       3,   # Branch;
                       5,   # DeletionN;
                       8,   # Nce;
                       10,  # Merge_A
                       11,  # Merge_C
                       12,  # Merge_G
                       13 ] # Merge_T

    fixedParamIdx = [ i for i in xrange(18) if i not in freeParamIdx ]
    fixedParamMask = [ (i in fixedParamIdx) for i in xrange(18) ]
    numFreeParams = len(freeParamIdx)

    fullStart = -0.23*np.array(~np.array(fixedParamMask), dtype=np.float32)

    """
    Starting point for training.
    """
    start = fullStart[freeParamIdx]
