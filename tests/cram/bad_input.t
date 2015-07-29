
Test for sane behavior in the presence of bogus arguments.

  $ variantCaller fake.alignmentset.xml -r fake.referenceset.xml -o test.fasta 2>&1 | tail -1
  variantCaller: error: Input file */fake.alignmentset.xml not found. (glob)
