
Check the gzip compressed output.

  $ export DATA=$TESTDIR/../data
  $ export INPUT=$DATA/lambda/aligned_reads_1.cmp.h5
  $ export REFERENCE=$DATA/lambda/lambdaNEB.fa
  $ variantCaller.py --algorithm=plurality -q 10 -r $REFERENCE $INPUT -o variants.gff.gz
  $ variantCaller.py --algorithm=plurality -q 10 -r $REFERENCE $INPUT -o consensus.csv.gz
  $ variantCaller.py --algorithm=plurality -q 10 -r $REFERENCE $INPUT -o consensus.fa.gz
  $ variantCaller.py --algorithm=plurality -q 10 -r $REFERENCE $INPUT -o consensus.fq.gz
  $ gunzip variants.gff.gz
  $ head variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.4
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus * (glob)
  ##source-commandline * (glob)
  ##sequence-region lambda_NEB3011 1 48502
  lambda_NEB3011\t.\tinsertion\t119\t119\t.\t.\t.\tvariantSeq=G;coverage=2;confidence=15;frequency=2;length=1 (esc)
  lambda_NEB3011\t.\tdeletion\t4517\t4517\t.\t.\t.\treference=T;coverage=2;confidence=15;frequency=2;length=1 (esc)
  lambda_NEB3011\t.\tinsertion\t6143\t6143\t.\t.\t.\tvariantSeq=T;coverage=2;confidence=15;frequency=2;length=1 (esc)
  $ gunzip consensus.csv.gz
  $ head consensus.csv
  referenceId,referencePos,coverage,consensus,consensusConfidence,consensusFrequency
  lambda_NEB3011,0,2,G,15,2
  lambda_NEB3011,1,2,G,15,2
  lambda_NEB3011,2,2,G,15,2
  lambda_NEB3011,3,2,C,15,2
  lambda_NEB3011,4,2,G,15,2
  lambda_NEB3011,5,2,G,3,1
  lambda_NEB3011,6,2,C,15,2
  lambda_NEB3011,7,2,GG,3,1
  lambda_NEB3011,8,2,A,15,2

