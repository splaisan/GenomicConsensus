
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/mruber/aligned_reads.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/mruber/Mruber_DSM_1279.fasta
  $ variantCaller.py -j${JOBS-8} --algorithm=quiver $INPUT -r $REFERENCE -o variants.gff -o css.fasta
