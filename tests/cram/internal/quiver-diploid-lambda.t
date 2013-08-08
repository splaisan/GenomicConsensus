
Reads are from a simulated diploid lambda, where there is a SNP at
each position 250 + 500k, and the SNP is a substitution "ACGT" ->
"CGTA".  How well do we pick up these SNPs?

  $ alias untabify="sed 's/\t/ /g'"
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/lambdaDiploid/aln.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/lambdaDiploid/lambda.fasta
  $ export EXPECTED_VARIANTS=/mnt/secondary/Share/Quiver/TestData/lambdaDiploid/v-expected.gff

  $ quiver --diploid $INPUT -r $REFERENCE \
  > -o variants.gff -o css.fasta

Take a look at the variants

  $ head -10 variants.gff | grep -v "#" | untabify
  lambda_NEB3011 . substitution 250 250 . . . reference=A;variantSeq=A/C;coverage=100;confidence=40
  lambda_NEB3011 . substitution 750 750 . . . reference=T;variantSeq=A/T;coverage=100;confidence=40
  lambda_NEB3011 . substitution 1250 1250 . . . reference=G;variantSeq=G/T;coverage=100;confidence=40

Use gffsubtract.pl to compare variants to expected.  Note that the
gffsubtract tool just looks at the coordinates, not the actual content
of the event, so it's not going to see if we called G/C as G/T, for
example.  Would be good to either write a better tool or make an easy
way to do this in Python.


False negatives:

  $ gffsubtract.pl $EXPECTED_VARIANTS variants.gff | untabify


False positives:

  $ gffsubtract.pl variants.gff $EXPECTED_VARIANTS | untabify
