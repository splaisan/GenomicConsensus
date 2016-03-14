Bite-sized quiver test using an All4Mers template!

  $ export DATA=$TESTDIR/../data
  $ export INPUT=$DATA/all4mer/out.aligned_subreads.bam
  $ export REFERENCE=$DATA/all4mer/All4mer.V2.01_Insert.fa

Run quiver.

  $ quiver $INPUT -r $REFERENCE -o v.gff -o css.fa


No variants!

  $ egrep -v '^#' v.gff | cat

Perfect consensus, no no-calls.

  $ cat css.fa
  >All4mer.V2.01_Insert|quiver
  CATCAGGTAAGAAAGTACGATGCTACAGCTTGTGACTGGTGCGGCACTTTTGGCTGAGTT
  TCCTGTCCACCTCATGTATTCTGCCCTAACGTCGGTCTTCACGCCATTACTAGACCGACA
  AAATGGAAGCCGGGGCCTTAAACCCCGTTCGAGGCGTAGCAAGGAGATAGGGTTATGAAC
  TCTCCCAGTCAATATACCAACACATCGTGGGACGGATTGCAGAGCGAATCTATCCGCGCT
  CGCATAATTTAGTGTTGATC
