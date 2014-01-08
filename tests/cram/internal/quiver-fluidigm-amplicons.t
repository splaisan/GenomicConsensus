
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/fluidigmAmplicons/aligned_reads.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/fluidigmAmplicons/MET_EGFR_Full_Genes.fasta
  $ quiver --noEvidenceConsensusCall=nocall -j${JOBS-8} $INPUT -r $REFERENCE \
  > -o variants.gff -o css.fasta
