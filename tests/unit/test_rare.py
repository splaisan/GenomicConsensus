import sys
from collections import defaultdict, OrderedDict
from nose import with_setup
from nose.tools import assert_equal, assert_true
from GenomicConsensus.options import options, parseOptions, importAdditionalDefaultOptions
from GenomicConsensus.rare import *
from GenomicConsensus.plurality.plurality import PluralityLocusSummary
from AlignmentHitStubs import *

class TestRareVariants(object):
    """Battery of tests for the rare variant algorithm. """

    snpsfreqs = ['T',0.0095,   # 19, < 1%
                 'G',0.025,    # 50, lower QV rare variant
                 '-',0.05,     # 100, highest QV rare variant
                 'A',0.9155]   # 1831, dominant allele (putative wild-type)

    def setup(self):
        sys.argv = ['','']
        parseOptions(relax=True)
        importAdditionalDefaultOptions(rare.additionalDefaultOptions)

    def test_consensus(self):
        """[RareAlignmentColumn.consensus] The core rare variant algo. """
        algCol = RareAlignmentColumn('test',1,'A')
        sf = TestRareVariants.snpsfreqs
        for i in xrange(0,8,2):
            count = int(sf[i+1]*2000)
            [algCol.addReadSnippet(sf[i]) for x in xrange(0, count)]

        results = algCol.consensus()

        # running just consensus will find any and all variants above 1% and
        # coverage > 500
        assert_equal(3, len(results))
        assert_equal(results[0], PluralityLocusSummary('test',1,2000,'G',84,50))
        assert_equal(results[1], PluralityLocusSummary('test',1,2000,'-',93,100))
        assert_equal(results[2], PluralityLocusSummary('test',1,2000,'A',93,1831))

    def gen_hits(self, coverage=2000):
        ref = "ATGG CTCC GATC TCAG".replace(' ','')
        hits = []
        sf = TestRareVariants.snpsfreqs
        for i in xrange(0,8,2):
            read = ref[0:9] + sf[i] + ref[10:]
            count = int(sf[i+1]*coverage)
            hit = AlignmentHitStub(0, FORWARD, ref, read)
            [hits.append(hit) for x in xrange(0, count)]

        return hits

    def test_coverage_option(self):
        """Can we adjust the minimum coverage"""
        options.variantCoverageThreshold = 110
        results = RareCaller.rare(('test',0,20), self.gen_hits(coverage=120))
        assert_equal(3, len(results))
        expected = [(('test',9),PluralityLocusSummary('test',9,119,'G',15,3)),
                    (('test',9),PluralityLocusSummary('test',9,119,'',37,6)),
                    (('test',9),PluralityLocusSummary('test',9,119,'A',93,109))]

        expected = set(expected)

        for r in results:
            assert_true(r in expected, "Unexpected result: %s" % str(r))


    def test_rare(self):
        """[RareCaller.rare] The rare variant alignment processing. """
        results = RareCaller.rare(('test',0,20),self.gen_hits())
        assert_equal(3, len(results))
        results.index((('test',9),PluralityLocusSummary('test',9,2000,'G',84,50)))
        results.index((('test',9),PluralityLocusSummary('test',9,2000,'',93,100)))
        results.index((('test',9),PluralityLocusSummary('test',9,2000,'A',93,1831)))
