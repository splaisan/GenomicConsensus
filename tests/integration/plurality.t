
Run plurality on the small example file, and make sure the GFF and
CSV output is correct.

  $ export DATA=$TESTDIR/../data   
  $ export INPUT=$DATA/lambda/aligned_reads_1.cmp.h5
  $ export REFERENCE=$DATA/lambda/lambdaNEB.fa
  $ variantCaller.py --algorithm=plurality -q 10 -r $REFERENCE -o variants.gff -o consensus.csv $INPUT

I like to show the head of the output files inline here so that glaringly obvious changes will
pop right out, but I verify that the files are exactly correct by looking at the md5 sums.

First, the variants.gff:

  $ head variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.3.3
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus v0.2.0
  ##source-commandline * (glob)
  ##sequence-header ref000001 lambda_NEB3011
  ##sequence-region ref000001 1 48502
  ref000001\t.\tinsertion\t119\t119\t.\t.\t.\tvariantSeq=G;coverage=2;confidence=15;frequency=2;length=1 (esc)
  ref000001\t.\tdeletion\t4517\t4517\t.\t.\t.\treference=T;coverage=2;confidence=15;frequency=2;length=1 (esc)

We expect 17 insertions, 8 deletions, and 0 SNVs.

  $ grep insertion variants.gff | wc | awk '{print $1}'
  17
  $ grep deletion variants.gff | wc | awk '{print $1}'
  8
  $ grep SNV variants.gff | wc | awk '{print $1}'
  0

Now make sure nothing else has changed on us!

  $ grep -v '\#.*' variants.gff | md5sum
  4322c060d417e904abf1beba21a612c3  -


Now, the consensus.csv:

  $ head consensus.csv
  referenceId,referencePos,coverage,consensus,consensusConfidence,consensusFrequency
  ref000001,0,2,G,15,2
  ref000001,1,2,G,15,2
  ref000001,2,2,G,15,2
  ref000001,3,2,C,15,2
  ref000001,4,2,G,15,2
  ref000001,5,2,G,3,1
  ref000001,6,2,C,15,2
  ref000001,7,2,GG,3,1
  ref000001,8,2,A,15,2
  $ md5sum consensus.csv
  ecdf6213d7132b5c178854ba35ad69c5  consensus.csv

