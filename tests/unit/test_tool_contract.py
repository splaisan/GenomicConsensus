
# TODO add quiver test (replacing cram test)

import unittest
import os.path

import pbcommand.testkit

DATA_DIR = os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")
assert os.path.isdir(DATA_DIR)

class TestGffToBed(pbcommand.testkit.PbTestApp):
    DRIVER_BASE = "gffToBed "
    DRIVER_EMIT = DRIVER_BASE + " --emit-tool-contract "
    DRIVER_RESOLVE = DRIVER_BASE + " --resolved-tool-contract "
    REQUIRES_PBCORE = True
    INPUT_FILES = [
        os.path.join(DATA_DIR, "converters", "variants.gff.gz"),
    ]
    TASK_OPTIONS = {
        "genomic_consensus.task_options.gff2bed_purpose": "variants",
        "genomic_consensus.task_options.track_name": None,
        "genomic_consensus.task_options.track_description": None,
        "genomic_consensus.task_options.use_score": 0,
    }

if __name__ == "__main__":
    unittest.main()
