
  $ export INPUT=/mnt/secondary/Quiver/TestData/lambda/job_038537.cmp.h5
  $ export REFERENCE=/mnt/secondary/Quiver/TestData/lambda/lambda.fasta
  $ variantCaller.py -j4 --algorithm=quiver $INPUT -r $REFERENCE -o variants.gff -o css.fasta
  $ cat variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.4
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus v0.2.0
  ##source-commandline * (glob)
  ##sequence-header ref000001 ref000001|lambda_NEB3011
  ##sequence-region ref000001 1 48502
