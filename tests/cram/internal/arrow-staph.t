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

Arrow, since the SNR capping fix, also gets no errors. 

  $ gffsubtract.pl arrow-variants.gff $MASK | grep -v "#" | sed 's/\t/ /g'

  $ fastacomposition arrow-css.fasta
  arrow-css.fasta A 960233 C 470725 G 470270 T 971457
