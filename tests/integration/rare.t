(If scipy isn't installed, we won't be able to run 'rare', so we should skip this test by
exiting 80))

  $ python -c "import scipy" > /dev/null 2>&1 || exit 80

Run plurality on the small example file, and make sure the GFF and
CSV output is correct.

  $ export DATA=$TESTDIR/../data   
  $ export INPUT=$DATA/bcrabl/aligned_reads_dccs.cmp.h5
  $ export REFERENCE=$DATA/bcrabl/bcrabl_amplicon_ucsf.fa
  $ variantCaller.py --algorithm=rare -r $REFERENCE -o variants.gff $INPUT

I like to show the head of the output files inline here so that glaringly obvious changes will
pop right out, but I verify that the files are exactly correct by looking at the md5 sums.

First, the variants.gff:

  $ cat variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.4
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus * (glob)
  ##source-commandline * (glob)
  ##sequence-region bcr_abl_amplicon 1 862
  bcr_abl_amplicon\t.\tsubstitution\t99\t99\t.\t.\t.\tvariantSeq=C;reference=T;coverage=2964;confidence=93;frequency=197;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t285\t285\t.\t.\t.\tvariantSeq=T;reference=A;coverage=2971;confidence=93;frequency=193;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t285\t285\t.\t.\t.\tvariantSeq=G;reference=A;coverage=2971;confidence=93;frequency=1307;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t286\t286\t.\t.\t.\tvariantSeq=G;reference=C;coverage=2971;confidence=93;frequency=125;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t286\t286\t.\t.\t.\tvariantSeq=T;reference=C;coverage=2971;confidence=93;frequency=716;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t297\t297\t.\t.\t.\tvariantSeq=G;reference=A;coverage=2971;confidence=93;frequency=211;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t301\t301\t.\t.\t.\tvariantSeq=G;reference=A;coverage=2971;confidence=4;frequency=30;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t406\t406\t.\t.\t.\tvariantSeq=G;reference=A;coverage=2975;confidence=93;frequency=309;length=1 (esc)
  bcr_abl_amplicon\t.\tsubstitution\t408\t408\t.\t.\t.\tvariantSeq=G;reference=A;coverage=2975;confidence=5;frequency=32;length=1 (esc)

