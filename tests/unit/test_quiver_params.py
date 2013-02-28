#
# Test aspects of the loading of quiver parameter sets from .ini files
#

from numpy.testing import assert_array_almost_equal
from nose.tools import assert_equal
import operator, numpy as np
from os.path import dirname as up

from GenomicConsensus.quiver.model import *


class StubParameterSet(object):
    def __init__(self, name):
        chem, modelName = name.split(".")[:2]
        if    modelName=="AllQVsModel":    model = AllQVsModel
        elif  modelName=="NoMergeQVModel": model = NoMergeQVModel
        elif  modelName=="NoQVsModel":     model = NoQVsModel
        self.name         = name
        self.chemistry    = chem
        self.model        = model
        self.quiverConfig = None


class TestBestParameterSet:
    def setup(self):
        self.parameterSets = { name : StubParameterSet(name)
                               for name in ["C5.AllQVsModel",
                                            "C5.NoMergeQVModel",
                                            "C5.NoQVsModel",
                                            "C6.NoQVsModel",
                                            "unknown.AllQVsModel",
                                            "unknown.NoMergeQVModel",
                                            "unknown.NoQVsModel"] }

    def test_bestParameterSet_1(self):
        assert_equal("C5.NoQVsModel",
                     bestParameterSet(self.parameterSets.values(), "C5", []).name)
        assert_equal("C5.NoMergeQVModel",
                     bestParameterSet(self.parameterSets.values(), "C5",
                                      NoMergeQVModel.requiredFeatures).name)
        assert_equal("C5.AllQVsModel",
                     bestParameterSet(self.parameterSets.values(), "C5",
                                      AllQVsModel.requiredFeatures).name)

    def test_bestParameterSet_2(self):
        # Try C6, where there is only one trained parameter set
        assert_equal("C6.NoQVsModel",
                     bestParameterSet(self.parameterSets.values(), "C6", []).name)
        assert_equal("C6.NoQVsModel",
                     bestParameterSet(self.parameterSets.values(), "C6",
                                      NoMergeQVModel.requiredFeatures).name)
        assert_equal("C6.NoQVsModel",
                     bestParameterSet(self.parameterSets.values(), "C6",
                                      AllQVsModel.requiredFeatures).name)

    def test_bestParameterSet_3(self):
        # Try the "future" chemistry C7, where we don't have access to
        # any trained parameter sets yet
        assert_equal("unknown.NoQVsModel",
                     bestParameterSet(self.parameterSets.values(), "C7", []).name)
        assert_equal("unknown.NoMergeQVModel",
                     bestParameterSet(self.parameterSets.values(), "C7",
                                      NoMergeQVModel.requiredFeatures).name)
        assert_equal("unknown.AllQVsModel",
                     bestParameterSet(self.parameterSets.values(), "C7",
                                      AllQVsModel.requiredFeatures).name)


class TestLoadingBundledParameters:

    def test_loadBundledParameterSet(self):
        """
        Make sure that the bundled parameter set in resources/ works.
        """
        paramSets = loadParameterSets(findParametersFile())
        assert "C2.AllQVsModel" in paramSets


    def test_loadParameterSetFromDir1(self):
        """
        Try loading from X, where the FS looks like:
          X/
            2013-03/
              GenomicConsensus/
                QuiverParameters.ini
        """
        quiverParamsIni = findParametersFile()
        X = up(up(up(quiverParamsIni)))
        paramSets = loadParameterSets(findParametersFile(X))
        assert "C2.AllQVsModel" in paramSets

    def test_loadParameterSetFromDir2(self):
        """
        Try loading from specified bundle
          X/
            2013-03/  <------------ here
              GenomicConsensus/
                QuiverParameters.ini
        """
        quiverParamsIni = findParametersFile()
        X = up(up(quiverParamsIni))
        paramSets = loadParameterSets(findParametersFile(X))
        assert "C2.AllQVsModel" in paramSets
