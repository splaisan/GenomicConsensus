
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/fluidigmAmplicons/aligned_reads.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/fluidigmAmplicons/MET_EGFR_Full_Genes.fasta
  $ quiver --noEvidenceConsensusCall=nocall -j${JOBS-8} $INPUT -r $REFERENCE \
  > -o variants.gff -o css.fasta
  [WARNING] This .cmp.h5 file lacks some of the QV data tracks that are required for optimal performance of the Quiver algorithm.  For optimal results use the ResequencingQVs workflow in SMRTPortal with bas.h5 files from an instrument using software version 1.3.1 or later, or the --forQuiver option to pbalign.
