
Check the gzip compressed output.

  $ export PATH=$TESTDIR/..:$PATH
  $ export INPUT=$TESTDIR/data/lambda/aligned_reads_1.cmp.h5
  $ export REFERENCE=$TESTDIR/data/lambda/lambdaNEB.fa
  $ variantCaller.py --algorithm=plurality -t 10 -r $REFERENCE $INPUT -o variants.gff.gz
  $ variantCaller.py --algorithm=plurality -t 10 -r $REFERENCE $INPUT -o consensus.csv.gz
  $ variantCaller.py --algorithm=plurality -t 10 -r $REFERENCE $INPUT -o consensus.fa.gz
  $ variantCaller.py --algorithm=plurality -t 10 -r $REFERENCE $INPUT -o consensus.fq.gz
  $ gunzip variants.gff.gz
  $ head variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.3.3
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus v0.2.0
  ##source-commandline * (glob)
  ##sequence-header ref000001 lambda_NEB3011
  ##sequence-region ref000001 1 48502
  ref000001\t.\tinsertion\t119\t119\t.\t.\t.\tlength=1;variantSeq=G;frequency=2;confidence=15;coverage=2 (esc)
  ref000001\t.\tinsertion\t1101\t1101\t.\t.\t.\tlength=1;variantSeq=C;frequency=2;confidence=15;coverage=2 (esc)
  $ gunzip consensus.csv.gz
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

