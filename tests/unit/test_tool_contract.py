
import unittest
import os.path

from pbcore.io import openDataSet, ContigSet
import pbcore.data
import pbcommand.testkit

# XXX local data directory, absolutely required
DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
assert os.path.isdir(DATA_DIR)

# optional (but required for TestSummarizeConsensus)
DATA_DIR_2 = "/mnt/secondary/Share/Quiver/TestData/tinyLambda/"

class TestQuiver(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "variantCaller "
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        pbcore.data.getBamAndCmpH5()[0],
#        "/pbi/dept/secondary/siv/testdata/SA3-DS/lambda/2372215/0007_tiny/Alignment_Results/m150404_101626_42267_c100807920800000001823174110291514_s1_p0.1.alignmentset.xml",
        pbcore.data.getLambdaFasta()
    ]
    TASK_OPTIONS = {
      "genomic_consensus.task_options.min_coverage": 0,
      "genomic_consensus.task_options.min_confidence": 0,
      "genomic_consensus.task_options.algorithm": "quiver",
      "genomic_consensus.task_options.diploid": False,
    }

    def run_after(self, rtc, output_dir):
        contigs_file = rtc.task.output_files[1]
        with openDataSet(contigs_file, strict=True) as ds:
            self.assertTrue(isinstance(ds, ContigSet))


class TestGffToBed(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "gffToBed "
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        os.path.join(DATA_DIR, "converters", "variants.gff.gz"),
    ]
    TASK_OPTIONS = {
        "genomic_consensus.task_options.track_name": "None",
        "genomic_consensus.task_options.track_description": "None",
        "genomic_consensus.task_options.use_score": 0,
    }


class TestGffToVcf(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "gffToVcf"
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        os.path.join(DATA_DIR, "converters", "variants.gff.gz"),
    ]
    TASK_OPTIONS = {
        "genomic_consensus.task_options.global_reference": "Staphylococcus_aureus_USA300_TCH1516",
    }


@unittest.skipUnless(os.path.isdir(DATA_DIR_2), "Missing %s" % DATA_DIR_2)
class TestSummarizeConsensus(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "summarizeConsensus"
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        os.path.join(DATA_DIR_2, "alignment_summary.gff"),
        os.path.join(DATA_DIR_2, "variants.gff.gz")
    ]
    TASK_OPTIONS = {}


if __name__ == "__main__":
    unittest.main()
