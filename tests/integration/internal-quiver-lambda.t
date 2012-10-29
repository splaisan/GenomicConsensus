
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/lambda/job_038537.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/lambda/lambda.fasta
  $ variantCaller.py -j8 --algorithm=quiver $INPUT -r $REFERENCE -o variants.gff -o css.fasta
  $ cat variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.4
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus * (glob)
  ##source-commandline * (glob)
  ##sequence-region lambda_NEB3011 1 48502
