Compare quiver vs. arrow on a high SNR Staph job.

  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/staphHighSnr/aligned_subreads.fofn
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/staph/S_aureus_USA300_TCH1516.fasta
  $ export MASK=/mnt/secondary/Share/Quiver/GenomeMasks/S_aureus_USA300_TCH1516-mask.gff

  $ quiver -j${JOBS-8} $INPUT -r $REFERENCE -o quiver-variants.gff -o quiver-css.fasta
  $ arrow  -j${JOBS-8} $INPUT -r $REFERENCE -o  arrow-variants.gff -o arrow-css.fasta

Quiver does a good job here---no errors.

  $ gffsubtract.pl quiver-variants.gff $MASK | grep -v "#" | sed 's/\t/ /g'

  $ fastacomposition quiver-css.fasta
  quiver-css.fasta A 960233 C 470725 G 470271 T 971458

Arrow, on the other hand, makes a number of errors.

  $ gffsubtract.pl arrow-variants.gff $MASK | grep -v "#" | sed 's/\t/ /g'
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 42450 42450 . . . reference=A;variantSeq=.;coverage=100;confidence=63
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 442861 442861 . . . reference=T;variantSeq=.;coverage=100;confidence=70
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 562107 562107 . . . reference=T;variantSeq=.;coverage=100;confidence=52
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 1041877 1041877 . . . reference=A;variantSeq=.;coverage=100;confidence=45
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 1801427 1801427 . . . reference=T;variantSeq=.;coverage=100;confidence=51
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 1840635 1840635 . . . reference=A;variantSeq=.;coverage=100;confidence=88
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 1850379 1850379 . . . reference=T;variantSeq=.;coverage=100;confidence=40
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 2147449 2147449 . . . reference=T;variantSeq=.;coverage=100;confidence=93
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 2231043 2231043 . . . reference=A;variantSeq=.;coverage=100;confidence=48
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 2587865 2587865 . . . reference=A;variantSeq=.;coverage=100;confidence=41

  $ fastacomposition arrow-css.fasta
  arrow-css.fasta A 960036 C 470641 G 470190 T 971288 a 178 c 84 g 80 t 158
