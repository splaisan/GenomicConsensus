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

import numpy as np, ConfigParser, collections, logging
from glob import glob
from os.path import join
from pkg_resources import resource_filename, Requirement

from GenomicConsensus.utils import die
from GenomicConsensus.quiver.utils import asFloatFeature, fst, snd
import ConsensusCore as cc

__all__ = [ "ParameterSet",
            "AllQVsModel",
            "NoMergeQVModel",
            "NoQVsModel",
            "AllQVsMergingByChannelModel",
            "NoQVsMergingByChannelModel",
            "QuiverConfig",
            "allQVsLoaded",
            "loadParameterSet",
            "loadQuiverConfig" ]

_basicParameterNames = \
          [ "Match"           ,
            "Mismatch"        , "MismatchS",
            "Branch"          , "BranchS",
            "DeletionN"       ,
            "DeletionWithTag" , "DeletionWithTagS",
            "Nce"             , "NceS",
            "Merge"           , "MergeS" ]

_mergeByChannelParameterNames = \
          [ "Match"           ,
            "Mismatch"        , "MismatchS",
            "Branch"          , "BranchS",
            "DeletionN"       ,
            "DeletionWithTag" , "DeletionWithTagS",
            "Nce"             , "NceS",
            "Merge_A"         , "Merge_C",
            "Merge_G"         , "Merge_T",
            "MergeS_A"        , "MergeS_C",
            "MergeS_G"        , "MergeS_T" ]

#
# Model classes
#

class Model(object):

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
        readIdentifier = str(aln.rowNumber)
        return cc.MappedRead(cls.extractFeatures(aln),
                             int(aln.RCRefStrand),
                             int(aln.referenceStart - windowStart),
                             int(aln.referenceEnd   - windowStart),
                             readIdentifier)

class AllQVsModel(Model):
    name = "AllQVsModel"
    rank = 3
    requiredFeatures = { "InsertionQV", "SubstitutionQV",
                         "DeletionQV" , "DeletionTag"   , "MergeQV" }
    parameterNames = _basicParameterNames

class NoMergeQVModel(Model):
    name = "NoMergeQVModel"
    rank = 2
    requiredFeatures = { "InsertionQV", "SubstitutionQV",
                         "DeletionQV" , "DeletionTag"   }
    parameterNames = _basicParameterNames

class NoQVsModel(Model):
    name = "NoQVsModel"
    rank = 1
    requiredFeatures = set([])
    parameterNames = _basicParameterNames

class AllQVsMergingByChannelModel(Model):
    name = "AllQVsMergingByChannelModel"
    rank = 4
    requiredFeatures = AllQVsModel.requiredFeatures
    parameterNames = _mergeByChannelParameterNames

class NoQVsMergingByChannelModel(Model):
    name = "NoQVsMergingByChannelModel"
    rank = -1
    requiredFeatures = set([])
    parameterNames = _mergeByChannelParameterNames


#
#   Code for accessing the ConsensusCore quiver parameter sets
#   from the .ini config file.
#

class ParameterSet(object):
    def __init__(self, name, model, chemistry, quiverConfig):
        self.name         = name
        self.chemistry    = chemistry
        self.model        = model
        self.quiverConfig = quiverConfig

def _getResourcesDirectory():
    return resource_filename(Requirement.parse("GenomicConsensus"),
                             "GenomicConsensus/quiver/resources")

def _majorityChemistry(cmpH5):
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

def _findParametersFile(filenameOrDirectory=None):
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

    #
    # Dirty hack for --diploid support and new scaling that I am using.
    # TODO: add banding options into quiver parameters file.
    #
    if (parameterSetName == "unknown.NoQVsModel" or
        parameterSetName == "P5-C3.AllQVsMergingByChannelModel"):
        bandingOptions     = cc.BandingOptions(4, 20)
        fastScoreThreshold = -50
    else:
        bandingOptions     = cc.BandingOptions(4, 5)
        fastScoreThreshold = -12.5

    quiverConfig = cc.QuiverConfig(qvModelParams,
                                   cc.ALL_MOVES,
                                   bandingOptions,
                                   fastScoreThreshold)
    return ParameterSet(parameterSetName, model, chem, quiverConfig)

def _loadParameterSets(iniFilename):
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

def _bestParameterSet(parameterSets, chemistry, qvsAvailable):
    fallbackParameterSets = \
        [ paramSet for paramSet in parameterSets.itervalues()
          if paramSet.chemistry == "unknown"
          if paramSet.model.requiredFeatures.issubset(qvsAvailable) ]
    perChemistryParameterSets = \
        [ paramSet for paramSet in parameterSets.itervalues()
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


#
#  QuiverConfig: the kitchen sink class of quiver options
#

class QuiverConfig(object):
    """
    Quiver configuration options
    """
    def __init__(self,
                 parameters,
                 minMapQV=10,
                 minPoaCoverage=3,
                 maxPoaCoverage=11,
                 mutationSeparation=10,
                 mutationNeighborhood=20,
                 maxIterations=40,
                 refineDinucleotideRepeats=True,
                 noEvidenceConsensus="nocall",
                 computeConfidence=True,
                 readStumpinessThreshold=0.1):

        self.minMapQV                   = minMapQV
        self.minPoaCoverage             = minPoaCoverage
        self.maxPoaCoverage             = maxPoaCoverage
        self.mutationSeparation         = mutationSeparation
        self.mutationNeighborhood       = mutationNeighborhood
        self.maxIterations              = maxIterations
        self.refineDinucleotideRepeats  = refineDinucleotideRepeats
        self.noEvidenceConsensus        = noEvidenceConsensus
        self.computeConfidence          = computeConfidence
        self.readStumpinessThreshold    = readStumpinessThreshold
        self.parameters                 = parameters

        # Convenience
        self.model                      = self.parameters.model
        self.ccQuiverConfig             = self.parameters.quiverConfig

    @staticmethod
    def _defaultQuiverParameters():
        return loadQuiverConfig("unknown.NoQVsModel")


#
#   Convenience functions
#

def allQVsLoaded(cmpH5):
    """
    Does this cmp.h5 file have the complete set of QV features?
    """
    return AllQVsModel.isCompatibleWithCmpH5(cmpH5)

def loadParameterSet(parameterSetNameOrCmpH5, parametersFile=None):
    parametersFile = _findParametersFile(parametersFile)
    sets = _loadParameterSets(parametersFile)
    if isinstance(parameterSetNameOrCmpH5, str):
        try:
            params = sets[parameterSetNameOrCmpH5]
        except:
            die("Quiver: no available parameter set named %s" % \
                parameterSetNameOrCmpH5)
    else:
        chemistry = _majorityChemistry(parameterSetNameOrCmpH5)
        qvsAvailable = parameterSetNameOrCmpH5.pulseFeaturesAvailable()
        params = _bestParameterSet(sets, chemistry, qvsAvailable)
    return params

def loadQuiverConfig(parameterSetNameOrCmpH5, parametersFile=None, **quiverConfigOpts):
    params = loadParameterSet(parameterSetNameOrCmpH5, parametersFile)
    return QuiverConfig(parameters=params, **quiverConfigOpts)
