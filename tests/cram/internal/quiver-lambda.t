Small lambda phage job, should be no errors.

  $ export CMPH5=/mnt/secondary/Share/Quiver/TestData/lambda.P4-C2/082796.cmp.h5
  $ export BAM=/mnt/secondary/Share/Quiver/TestData/lambda.P4-C2/082796.bam
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/lambda.P4-C2/lambdaNEB.fa

First try the cmp.h5:

  $ mkdir CmpH5; cd CmpH5
  $ quiver -j${JOBS-16} --noEvidenceConsensusCall=nocall $CMPH5 -r $REFERENCE \
  > -o variants.gff -o css.fasta.gz -o css.fq.gz

  $ grep -v "##" variants.gff
  [1]

  $ gunzip css.fasta.gz
  $ fastadiff -c FALSE css.fasta $REFERENCE

  $ cd ..

#Now run on the BAM:
#
#  $ mkdir BAM; cd BAM
#  $ quiver -j${JOBS-16} --noEvidenceConsensusCall=nocall $BAM -r $REFERENCE \
#  > -o variants.gff -o css.fasta.gz -o css.fq.gz
#  [WARNING] 'fancyChunking' not yet available for BAM, disabling
#
#  $ grep -v "##" variants.gff
#  [1]
#
#  $ gunzip css.fasta.gz
#  $ fastadiff -c FALSE css.fasta $REFERENCE
