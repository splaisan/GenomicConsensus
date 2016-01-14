
Test performance of summarizeConsensus with large numbers of variants.  This
should not take more than 5 minutes or so unless an O(N^2) loop is used.

  $ TESTDATA="/pbi/dept/secondary/siv/testdata/pbreports-unittest/data/summarizeConsensus"
  $ VARIANTS=$TESTDATA/variants.gff
  $ SUMMARY=$TESTDATA/alignment_summary.gff
  $ summarizeConsensus --variantsGff $VARIANTS --output consensus.gff $SUMMARY
