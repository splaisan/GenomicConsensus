
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/sequel-6k-ecoli/all_chunks.alignmentset.xml
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/ecoli/ecoliK12_pbi_March2013.fasta
  $ export MASK=/mnt/secondary/Share/Quiver/GenomeMasks/ecoliK12_pbi_March2013-mask.gff

  $ arrow  -j${JOBS-16} $INPUT -r $REFERENCE -o  arrow-variants.gff -o arrow-css.fasta

We do all right on short-insert Sequel data since Flea---only one
error identified here:

  $ gffsubtract.pl arrow-variants.gff $MASK | grep -v "#" | sed 's/\t/ /g'
  ecoliK12_pbi_March2013 . insertion 3468072 3468072 . . . reference=.;variantSeq=A;coverage=100;confidence=93
