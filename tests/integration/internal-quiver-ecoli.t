
  $ export INPUT=/mnt/secondary/Quiver/TestData/ecoli/job_044601.cmp.h5
  $ export REFERENCE=/mnt/secondary/Quiver/TestData/ecoli/ecoli_mutated.fasta
  $ qrsh -V -pe smp 16 -q secondary \
  >   variantCaller.py -j16 --algorithm=quiver $INPUT -r $REFERENCE -o variants.gff -o css.fasta
