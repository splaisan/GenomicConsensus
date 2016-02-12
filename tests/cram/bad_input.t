
Test for sane behavior in the presence of bogus arguments.

  $ variantCaller fake.alignmentset.xml -r fake.referenceset.xml -o test.fasta 2>&1 | tail -1
  variantCaller: error: Input file */fake.alignmentset.xml not found. (glob)

Test that it doesn't crash when one BAM file in an otherwise valid AlignmentSet is empty.

  $ DATA=$TESTDIR/../data/sanity
  $ REF="`python -c 'import pbcore.data ; print pbcore.data.getLambdaFasta()'`"
  $ variantCaller --reference $REF -o contig.fasta $DATA/mixed.alignmentset.xml
  [WARNING] */empty.subreads.bam contains no reads! (glob)
  [WARNING] */empty.subreads.bam contains no reads! (glob)
