#
# Test aspects of the loading of quiver parameter sets from .ini files
#

from nose.tools import assert_equal
from os.path import dirname as up

import GenomicConsensus.quiver.model as m

class StubParameterSet(object):
    def __init__(self, name):
        chem, modelName = name.split(".")[:2]
        if    modelName=="AllQVsModel":    model = m.AllQVsModel
        elif  modelName=="NoMergeQVModel": model = m.NoMergeQVModel
        elif  modelName=="NoQVsModel":     model = m.NoQVsModel
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
                     m._bestParameterSet(self.parameterSets, "C5", []).name)
        assert_equal("C5.NoMergeQVModel",
                     m._bestParameterSet(self.parameterSets, "C5",
                                      m.NoMergeQVModel.requiredFeatures).name)
        assert_equal("C5.AllQVsModel",
                     m._bestParameterSet(self.parameterSets, "C5",
                                      m.AllQVsModel.requiredFeatures).name)

    def test_bestParameterSet_2(self):
        # Try C6, where there is only one trained parameter set
        assert_equal("C6.NoQVsModel",
                     m._bestParameterSet(self.parameterSets, "C6", []).name)
        assert_equal("C6.NoQVsModel",
                     m._bestParameterSet(self.parameterSets, "C6",
                                        m.NoMergeQVModel.requiredFeatures).name)
        assert_equal("C6.NoQVsModel",
                     m._bestParameterSet(self.parameterSets, "C6",
                                        m.AllQVsModel.requiredFeatures).name)

    def test_bestParameterSet_3(self):
        # Try the "future" chemistry C7, where we don't have access to
        # any trained parameter sets yet
        assert_equal("unknown.NoQVsModel",
                     m._bestParameterSet(self.parameterSets, "C7", []).name)
        assert_equal("unknown.NoMergeQVModel",
                     m._bestParameterSet(self.parameterSets, "C7",
                                        m.NoMergeQVModel.requiredFeatures).name)
        assert_equal("unknown.AllQVsModel",
                     m._bestParameterSet(self.parameterSets, "C7",
                                        m.AllQVsModel.requiredFeatures).name)


class TestLoadingBundledParameters:

    def test_loadBundledParameterSet(self):
        """
        Make sure that the bundled parameter set in resources/ works.
        """
        paramSets = m._loadParameterSets(m._findParametersFile())
        assert "C2.AllQVsModel" in paramSets


    def test_loadParameterSetFromDir1(self):
        """
        Try loading from X, where the FS looks like:
          X/
            2013-03/
              GenomicConsensus/
                QuiverParameters.ini
        """
        quiverParamsIni = m._findParametersFile()
        X = up(up(up(quiverParamsIni)))
        paramSets = m._loadParameterSets(m._findParametersFile(X))
        assert "C2.AllQVsModel" in paramSets

    def test_loadParameterSetFromDir2(self):
        """
        Try loading from specified bundle
          X/
            2013-03/  <------------ here
              GenomicConsensus/
                QuiverParameters.ini
        """
        quiverParamsIni = m._findParametersFile()
        X = up(up(quiverParamsIni))
        paramSets = m._loadParameterSets(m._findParametersFile(X))
        assert "C2.AllQVsModel" in paramSets
