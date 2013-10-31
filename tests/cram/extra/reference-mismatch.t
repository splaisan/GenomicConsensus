
Test a few scenarios where the reference FASTA disagrees slightly from
the contigs aligned against in the cmp.h5, and make sure things behave
sanely.

  $ export DATA=$TESTDIR/../../data
  $ export INPUT=$DATA/hcv/aligned_reads.cmp.h5
  $ export WRONG_REFERENCE=$DATA/fluidigm_amplicons/Fluidigm_human_amplicons.fasta
  $ export REFERENCE_SUBSET=$DATA/hcv/5primeEnd.fa
  $ export REFERENCE_NO_FAI=$DATA/hcv/3primeEnd.fa

No .fai file:

  $ quiver $INPUT -r $REFERENCE_NO_FAI -o variants.gff -o consensus.fastq
  Companion FASTA index (.fai) file not found or malformatted! Use 'samtools faidx' to generate FASTA index.
  [255]

Wrong reference:

  $ quiver $INPUT -r $WRONG_REFERENCE -o variants.gff -o consensus.fastq
  No reference groups in the FASTA file were aligned against.  Did you select the wrong reference FASTA file?
  [255]

Reference containing a subset of the reference that was aligned to:

  $ quiver $INPUT -r $REFERENCE_SUBSET -o variants.gff -o consensus.fastq
  * [WARNING] Some reference contigs aligned against are not found in the reference FASTA.  Will process only those contigs supported by the reference FASTA. (glob)
  * [WARNING] This .cmp.h5 file lacks some of the QV data tracks that are required for optimal performance of the Quiver algorithm.  For optimal results use the ResequencingQVs workflow in SMRTPortal with bas.h5 files from an instrument using software version 1.3.1 or later. (glob)
