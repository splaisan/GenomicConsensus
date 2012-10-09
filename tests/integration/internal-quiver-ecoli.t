
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/ecoli/job_044601.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/ecoli/ecoli_mutated.fasta
  $ variantCaller.py -j8 --algorithm=quiver $INPUT -r $REFERENCE -o variants.gff -o css.fasta
