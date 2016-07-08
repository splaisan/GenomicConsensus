
import unittest
import os.path

from pbcore.io import openDataSet, ContigSet
import pbcommand.testkit

import pbtestdata

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
assert os.path.isdir(DATA_DIR)

# These tests seem to cause some logging failure at shutdown;
# disabling them pending upstream fix.  See:
# https://bugzilla.nanofluidics.com/show_bug.cgi?id=33699

# class TestVariantCaller(pbcommand.testkit.PbTestApp):
#     DRIVER_BASE = "variantCaller "
#     DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
#     DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
#     REQUIRES_PBCORE = True
#     INPUT_FILES = [
#         pbtestdata.get_file("aligned-xml"), pbtestdata.get_file("lambdaNEB")
#     ]
#     TASK_OPTIONS = {
#       "genomic_consensus.task_options.min_coverage": 0,
#       "genomic_consensus.task_options.min_confidence": 0,
#       "genomic_consensus.task_options.algorithm": "quiver",
#       "genomic_consensus.task_options.diploid": False,
#     }

#     def test_run_e2e(self):
#         import ipdb; ipdb.set_trace()
#         super(TestVariantCaller, self).test_run_e2e()

#     def run_after(self, rtc, output_dir):
#         contigs_file = rtc.task.output_files[1]
#         with openDataSet(contigs_file, strict=True) as ds:
#             self.assertTrue(isinstance(ds, ContigSet))


# class TestVariantCallerArrow(TestVariantCaller):
#     TASK_OPTIONS = {
#       "genomic_consensus.task_options.algorithm": "arrow",
#     }


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


class TestSummarizeConsensus(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "summarizeConsensus"
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        pbtestdata.get_file("alignment-summary-gff"),
        pbtestdata.get_file("variants-gff")
    ]
    TASK_OPTIONS = {}


if __name__ == "__main__":
    unittest.main()
