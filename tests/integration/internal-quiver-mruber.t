
  $ export INPUT=/mnt/secondary/Quiver/TestData/mruber/aligned_reads.cmp.h5
  $ export REFERENCE=/mnt/secondary/Quiver/TestData/mruber/MRuber_CA_Assembly.fasta
  $ variantCaller.py -j16 --algorithm=quiver $INPUT -r $REFERENCE -o variants.gff -o css.fasta
