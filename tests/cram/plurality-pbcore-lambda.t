
Run plurality on the lambda file in pbcore.

  $ export CMPH5=`python -c "from pbcore import data as D; print D.getBamAndCmpH5()[1]"`
  $ export BAM=`python -c "from pbcore import data as D; print D.getBamAndCmpH5()[0]"`
  $ export REF=`python -c "from pbcore import data as D; print D.getLambdaFasta()"`

  $ plurality $CMPH5 -r $REF -o css-cmph5.fa -o v-cmph5.gff
  $ cat v-cmph5.gff | grep -v "#" | sed 's/	/ /g'
  lambda_NEB3011 . deletion 4945 4945 . . . reference=C;variantSeq=.;frequency=4;coverage=7;confidence=40

  $ fastacomposition css-cmph5.fa
  css-cmph5.fa A 960 C 994 G 1129 T 941 a 11376 c 10366 g 11689 t 11047

  $ plurality $BAM -r $REF -o css-bam.fa -o v-bam.gff
  [WARNING] 'fancyChunking' not yet available for BAM, disabling
  $ cat v-bam.gff | grep -v "#" | sed  's/	/ /g'
  lambda_NEB3011 . deletion 4945 4945 . . . reference=C;variantSeq=.;frequency=4;coverage=7;confidence=40

  $ fastacomposition css-bam.fa
  css-bam.fa A 960 C 994 G 1129 T 941 a 11376 c 10366 g 11689 t 11047

  $ fastadiff css-cmph5.fa css-bam.fa
