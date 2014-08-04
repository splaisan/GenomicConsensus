  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/lambda/job_038537.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/lambda/lambda.fasta

  $ quiver -p C2 -j${JOBS-8} --noEvidenceConsensusCall=nocall $INPUT -r $REFERENCE \
  > -o variants.gff -o css.fasta.gz -o css.fq.gz

  $ grep -v "##" variants.gff
  [1]

  $ gunzip css.fasta.gz
  $ fastadiff -c FALSE css.fasta $REFERENCE
