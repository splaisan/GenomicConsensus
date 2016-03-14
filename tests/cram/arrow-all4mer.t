Bite-sized quiver test using an All4Mers template!

  $ export DATA=$TESTDIR/../data
  $ export INPUT=$DATA/all4mer/out.aligned_subreads.bam
  $ export REFERENCE=$DATA/all4mer/All4mer.V2.01_Insert.fa

Run arrow.

  $ arrow $INPUT -r $REFERENCE -o v.gff -o css.fa

Perfect consensus!

  $ fastadiff -c FALSE css.fa $REFERENCE

No variants!

  $ egrep -v '^#' v.gff | cat

No no-calls.

  $ fastacomposition css.fa
  css.fa A 65 C 66 G 64 T 65
