
Some tests of a "fluidigm amplicons" dataset

  $ export DATA=$TESTDIR/../data   
  $ export INPUT=$DATA/fluidigm_amplicons/040500.cmp.h5
  $ export REFERENCE=$DATA/fluidigm_amplicons/Fluidigm_human_amplicons.fasta

Set the QV threshold to 10.

  $ variantCaller.py --algorithm=plurality -r $REFERENCE -q 10 -o variants.gff -o consensus.csv -o consensus.fastq $INPUT

There are two true SNVs (and one diploid SNV that we miss right now).

  $ grep insertion variants.gff | wc | awk '{print $1}'
  0
  $ grep deletion variants.gff | wc | awk '{print $1}'
  0
  $ grep substitution variants.gff
  EGFR_Exon_23\t.\tsubstitution\t48\t48\t.\t.\t.\tvariantSeq=C;reference=T;coverage=248;confidence=93;frequency=239;length=1 (esc)


