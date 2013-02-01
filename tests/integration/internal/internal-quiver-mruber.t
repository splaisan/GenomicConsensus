
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/mruber/aligned_reads.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/mruber/MRuber_CA_Assembly.fasta
  $ variantCaller.py -j8 --algorithm=quiver $INPUT -r $REFERENCE -o variants.gff -o css.fasta
