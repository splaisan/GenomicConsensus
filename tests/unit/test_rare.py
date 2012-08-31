from collections import defaultdict
from nose.tools import assert_equal
from GenomicConsensus.rare.rare import *
from AlignmentHitStubs import *

class TestRareVariants:
    """Battery of tests for the rare variant algorithm. """

    snpsfreqs = ['T',19,   # < 1%
                 'G',69,   # lower QV rare variant
                 '-',169,  # highest QV rare variant
                 'A',2000] # dominant allele (putative wild-type) 

    def test_consensus(self):
        """[RareAlignmentColum.consensus] The core rare variant algo. """
        algCol = RareAlignmentColumn('test',1,'A')
        j = 0
        for i in range(1,2000):
            algCol.addReadSnippet(TestRareVariants.snpsfreqs[j]) 
            if i % TestRareVariants.snpsfreqs[j+1] == 0: j = j + 2

        results = algCol.consensus()

        # running just consensus will find any and all variants above 1% and 
        # coverage > 500
        assert_equal(3, len(results))
        assert_equal(results[0], PluralityLocusSummary('test',1,1999,'G',84,50))
        assert_equal(results[1], PluralityLocusSummary('test',1,1999,'-',93,100))
        assert_equal(results[2], PluralityLocusSummary('test',1,1999,'A',93,1830))

    def gen_hits(self):
        ref = "ATGG CTCC GATC TCAG".replace(' ','')
        hits = []
        j = 0
        for i in range(1,2000):
            read = ref[0:9] + TestRareVariants.snpsfreqs[j] + ref[10:]
            hits.append(AlignmentHitStub(0, FORWARD, ref, read))
            if i % TestRareVariants.snpsfreqs[j+1] == 0: j = j + 2

        return hits

    def test_rare(self):
        """[RareCaller.rare] The rare variant alignment processing. """
        results = RareCaller.rare(('test',0,20),self.gen_hits())
        assert_equal(3, len(results))
        results.index((('test',9),PluralityLocusSummary('test',9,1999,'G',84,50)))
        results.index((('test',9),PluralityLocusSummary('test',9,1999,'',93,100)))
        results.index((('test',9),PluralityLocusSummary('test',9,1999,'A',93,1830)))

    def test_consumer_entry(self):
        """Tests the consumer entry point. """
        # some local imports that we'll be 'mocking'
        from GenomicConsensus import reference
        from GenomicConsensus.io.consumers import consumer
        from GenomicConsensus.options import options

        # monkey patches to enable unit testing
        options.referenceWindow = ''
        options.referenceChunkSize = 1

        # reference.numChunks override
        def testNumChunks(refId, referenceChunkSize, referenceWindow=None):
            return 1

        reference.numChunks = testNumChunks
        reference.byId['test'] = reference.ReferenceEntry('test','test',
                                                          'test','AAAA',16)

        # now we're getting somewhere, define consumer so we can inspect on
        # the way in.  This is the actual unit test.  Note the deletion should
        # be gone now and only substitutions remain.
        @consumer
        def testConsumer(arg,**kwargs):
            while True:
                (refId, tbl) = (yield)
                assert_equal(dict,type(tbl))
                assert_equal(2, sum([len(v) for v in tbl.values()]))
                assert_equal(tbl[9][0], PluralityLocusSummary('test',9,1999,'G',84,50))
                assert_equal(tbl[9][1], PluralityLocusSummary('test',9,1999,'A',93,1830))

        # now import the result class ... after we hack the reference module
        from GenomicConsensus.rare.rare import RareResult
        rr = RareResult()

        # some dummy values
        rr.consensusByRefId = OrderedDict()
        rr.chunksReceivedById = defaultdict(int)
        
        # setup our sneaky consumer
        rr.consumers = [testConsumer(1,fn='blah')]

        # now start the music
        rc = RareCaller()
        result = rc.onChunk(('test',0,20),self.gen_hits())
        rr.onResult(result)

