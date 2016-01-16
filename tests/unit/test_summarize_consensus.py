
"""
Test summarizeConsensus with synthetic inputs, which should include all
boundary conditions.
"""

import subprocess
import tempfile
import unittest
import os


VARIANTS = """\
chr1\t.\tsubstitution\t100\t100\t.\t.\t.\treference=C;variantSeq=A;coverage=38;confidence=48
chr1\t.\tsubstitution\t5001\t5001\t.\t.\t.\treference=C;variantSeq=A;coverage=38;confidence=48
chr1\t.\tsubstitution\t10000\t10000\t.\t.\t.\treference=A;variantSeq=G;coverage=38;confidence=48
chr1\t.\tinsertion\t12000\t12000\t.\t.\t.\treference=.;variantSeq=G;coverage=38;confidence=48
chr1\t.\tdeletion\t15001\t15001\t.\t.\t.\treference=T;variantSeq=.;coverage=38;confidence=49
chr2\t.\tdeletion\t20000\t20000\t.\t.\t.\treference=T;variantSeq=.;coverage=38;confidence=49
chr2\t.\tsubstitution\t23469\t23469\t.\t.\t.\treference=T;variantSeq=A;coverage=38;confidence=47
chr3\t.\tsubstitution\t3469\t3469\t.\t.\t.\treference=C;variantSeq=A;coverage=32;confidence=40"""

SUMMARY = """\
chr1\t.\tregion\t1\t5000\t0.00\t+\t.\tcov=4,23,28;cov2=20.162,5.851;gaps=0,0
chr1\t.\tregion\t5001\t10000\t0.00\t+\t.\tcov=18,24,29;cov2=24.009,2.316;gaps=0,0
chr1\t.\tregion\t10001\t15000\t0.00\t+\t.\tcov=17,25,29;cov2=24.182,2.596;gaps=0,0
chr1\t.\tregion\t15001\t20000\t0.00\t+\t.\tcov=20,34,49;cov2=33.714,7.143;gaps=0,0
chr2\t.\tregion\t1\t5000\t0.00\t+\t.\tcov=18,49,80;cov2=49.018,10.125;gaps=0,0
chr2\t.\tregion\t10000\t23469\t0.00\t+\t.\tcov=0,48,89;cov2=47.303,12.036;gaps=1,24
chr3\t.\tregion\t1\t7000\t0.00\t+\t.\tcov=0,48,89;cov2=47.303,12.036;gaps=1,24"""

EXPECTED = """\
##source GenomicConsensus 2.0.0
##pacbio-alignment-summary-version 0.6
##source-commandline this line will be skipped in the comparison
chr1\t.\tregion\t1\t5000\t0.00\t+\t.\tcov=4,23,28;cov2=20.162,5.851;gaps=0,0;cQv=20,20,20;del=0;ins=0;sub=1
chr1\t.\tregion\t5001\t10000\t0.00\t+\t.\tcov=18,24,29;cov2=24.009,2.316;gaps=0,0;cQv=20,20,20;del=0;ins=0;sub=2
chr1\t.\tregion\t10001\t15000\t0.00\t+\t.\tcov=17,25,29;cov2=24.182,2.596;gaps=0,0;cQv=20,20,20;del=0;ins=1;sub=0
chr1\t.\tregion\t15001\t20000\t0.00\t+\t.\tcov=20,34,49;cov2=33.714,7.143;gaps=0,0;cQv=20,20,20;del=1;ins=0;sub=0
chr2\t.\tregion\t1\t5000\t0.00\t+\t.\tcov=18,49,80;cov2=49.018,10.125;gaps=0,0;cQv=20,20,20;del=0;ins=0;sub=0
chr2\t.\tregion\t10000\t23469\t0.00\t+\t.\tcov=0,48,89;cov2=47.303,12.036;gaps=1,24;cQv=20,20,20;del=1;ins=0;sub=1
chr3\t.\tregion\t1\t7000\t0.00\t+\t.\tcov=0,48,89;cov2=47.303,12.036;gaps=1,24;cQv=20,20,20;del=0;ins=0;sub=1"""

class TestSummarizeConsensus(unittest.TestCase):

    def setUp(self):
        self.variants_gff = tempfile.NamedTemporaryFile(suffix=".gff").name
        self.summary_gff = tempfile.NamedTemporaryFile(suffix=".gff").name
        with open(self.variants_gff, "w") as gff:
            gff.write(VARIANTS)
        with open(self.summary_gff, "w") as gff:
            gff.write(SUMMARY)

    def tearDown(self):
        os.remove(self.variants_gff)
        os.remove(self.summary_gff)

    def test_integration(self):
        gff_out = tempfile.NamedTemporaryFile(suffix=".gff").name
        args = [
            "summarizeConsensus",
            "--variants", self.variants_gff,
            "--output", gff_out,
            self.summary_gff
        ]
        self.assertEqual(subprocess.call(args), 0)
        with open(gff_out) as gff:
            lines = gff.read().splitlines()
            expected_lines = EXPECTED.splitlines()  
            self.assertEqual(len(lines), len(expected_lines))
            for a, b in zip(lines, expected_lines):
                if not a.startswith("##source-commandline"):
                    self.assertEqual(a, b)


if __name__ == "__main__":
    unittest.main()
