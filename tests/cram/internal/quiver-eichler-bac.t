
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/eichler/053727.cmp.h5
  $ export SANGER_REFERENCE=/mnt/secondary/Share/Quiver/TestData/eichler/CH17-157L1.finished.fa
  $ export ASSEMBLY_REFERENCE=/mnt/secondary/Share/Quiver/TestData/eichler/CH17_157L1_quiver_fasta.fasta

The QVs warning gets printed to stderr N times ... ignore it for now.

  $ quiver --noEvidenceConsensusCall=nocall \
  > -j${JOBS-8} $INPUT -r $ASSEMBLY_REFERENCE -o variants.gff -o css.fasta 2>/dev/null

Variant scores are currently miscalibrated (need to fix the
NoMergeQVModel; bug 22255).  Note that these variants listed below are
reckoned compared to the assembly reference, so they are not really
variants so much as errors in the assembly.  Variants assessed using
MuMMer at the end are compared to the Sanger reference.

  $ sed 's/\t/ /g' variants.gff | grep -v '#'
  CH17-157L1 . deletion 797 797 . . . reference=G;variantSeq=.;coverage=100;confidence=48
  CH17-157L1 . deletion 805 805 . . . reference=T;variantSeq=.;coverage=100;confidence=47
  CH17-157L1 . deletion 26174 26175 . . . reference=AC;variantSeq=.;coverage=100;confidence=48
  CH17-157L1 . deletion 93356 93357 . . . reference=CG;variantSeq=.;coverage=100;confidence=48
  CH17-157L1 . insertion 230679 230679 . . . reference=.;variantSeq=A;coverage=100;confidence=49
  CH17-157L1 . insertion 230681 230681 . . . reference=.;variantSeq=CA;coverage=100;confidence=48
  CH17-157L1 . insertion 230684 230684 . . . reference=.;variantSeq=C;coverage=100;confidence=48


  $ fastacomposition css.fasta
  css.fasta A 65434 C 51200 G 50136 N 1028 T 63110

Use the MuMMer suite to look at the differences from the reference.

  $ nucmer -mum $SANGER_REFERENCE css.fasta 2>/dev/null

First: no structural differences.

  $ show-diff -H -q out.delta | sed 's/\t/ /g'
  CH17-157L1|quiver BRK 1 517 517
  CH17-157L1|quiver GAP 55509 56008 500 500 0
  CH17-157L1|quiver BRK 230889 230908 20

Next, the SNPs.

  $ show-snps -H -C -x10 out.delta
     16035   . A   16054     |     8523    16035  |  AAAAAAAAAA.ATTGCTTGCC  AAAAAAAAAAAATTGCTTGCC  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     24558   . A   24578     |     8523    24558  |  AAAAAAAAAA.AGCCTGGATG  AAAAAAAAAAAAGCCTGGATG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     51215   C .   51234     |     4275    51215  |  GGCCCGCCCCCCGGGCAGCCA  GGCCCGCCCC.CGGGCAGCCA  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     64634   C .   64652     |     8645    64634  |  GACCCCCCCCCCACCGGTCAG  GACCCCCCCC.CACCGGTCAG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     85478   . T   85497     |     8834    85478  |  TTTTTTTTTT.TACTAACCAG  TTTTTTTTTTTTACTAACCAG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     94312   . T   94332     |     8834    94312  |  TTTTTTTTTT.TAGACAGAGT  TTTTTTTTTTTTAGACAGAGT  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
    106985   . T   107006    |        0   106985  |  TTTTTTTTTT.TCCTGAGCAG  TTTTTTTTTTTTTCCTGAGCA  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
    106985   . T   107007    |        0   106985  |  TTTTTTTTTT.TCCTGAGCAG  TTTTTTTTTTTTCCTGAGCAG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
    182920   . A   182943    |    47946    47946  |  AAAAAAAAAA.ATGTGGTCTC  AAAAAAAAAAAATGTGGTCTC  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
