
Test performance of summarizeConsensus with large numbers of variants.

  $ TESTDATA="/pbi/dept/secondary/siv/testdata/pbreports-unittest/data/summarizeConsensus"

First test has 48000 regions total but hundreds of small contigs, with more
than 5 million variants.  This should not take more than 5 minutes or so
unless an O(N^2) loop is used.

  $ VARIANTS=$TESTDATA/variants.gff
  $ SUMMARY=$TESTDATA/alignment_summary.gff
  $ GFFREF=$TESTDATA/alignment_summary_variants.gff
  $ OUTPUT=alignment_summary_variants_test.gff
  $ summarizeConsensus --variantsGff $VARIANTS --output $OUTPUT $SUMMARY
  $ diff -I "##source" $OUTPUT $GFFREF

Second test has 20000 regions in a single contig, and 10000 variants.  This
will also take several minutes.

  $ VARIANTS=$TESTDATA/variants_big_chr.gff
  $ SUMMARY=$TESTDATA/alignment_summary_big_chr.gff
  $ GFFREF=$TESTDATA/alignment_summary_variants_big_chr.gff
  $ OUTPUT=alignment_summary_variants_big_chr_test.gff
  $ summarizeConsensus --variantsGff $VARIANTS --output $OUTPUT $SUMMARY
  $ diff -I "##source" $OUTPUT $GFFREF
