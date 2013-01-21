
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

  $ sed 's/\t/ /g' variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.4
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus * (glob)
  ##source-commandline * (glob)
  ##sequence-region CH17-157L1 1 230921
  CH17-157L1 . deletion 797 797 . . . reference=G;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 805 805 . . . reference=T;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 10965 10965 . . . reference=T;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 15160 15160 . . . reference=A;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 24556 24556 . . . reference=A;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 26174 26175 . . . reference=AC;coverage=100;confidence=93;length=2
  CH17-157L1 . deletion 51238 51238 . . . reference=C;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 64654 64654 . . . reference=C;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 93356 93357 . . . reference=CG;coverage=100;confidence=93;length=2
  CH17-157L1 . deletion 161413 161413 . . . reference=T;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 180807 180807 . . . reference=T;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 190890 190890 . . . reference=T;coverage=100;confidence=93;length=1
  CH17-157L1 . deletion 202304 202304 . . . reference=T;coverage=100;confidence=93;length=1


  $ fastacomposition css.fasta
  css.fasta A 65472 C 51211 G 50142 N 921 T 63159

Use the MuMMer suite to look at the differences from the reference.

  $ nucmer -mum $SANGER_REFERENCE css.fasta 2>/dev/null

First: no structural differences.

  $ show-diff -q out.delta | sed 's/\t/ /g'
  * (glob)
  NUCMER
  
  [SEQ] [TYPE] [S1] [E1] [LEN 1]
  CH17-157L1|quiver BRK 1 500 500
  CH17-157L1|quiver BRK 230485 230905 421

Next, the SNPs.

  $ show-snps -H -C -x10 out.delta
     24558   . A   24577     |    24076    24558  |  AAAAAAAAAA.AGCCTGGATG  AAAAAAAAAAAAGCCTGGATG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     51215   C .   51233     |     1765    51215  |  GGCCCGCCCCCCGGGCAGCCA  GGCCCGCCCC.CGGGCAGCCA  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     52980   . A   52999     |     1765    52980  |  AAAAAAAAAA.ACAACAAACA  AAAAAAAAAAAACAACAAACA  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     64634   C .   64652     |    11654    64634  |  GACCCCCCCCCCACCGGTCAG  GACCCCCCCC.CACCGGTCAG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     85478   . T   85497     |     8834    85478  |  TTTTTTTTTT.TACTAACCAG  TTTTTTTTTTTTACTAACCAG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
     94312   . T   94332     |     8834    94312  |  TTTTTTTTTT.TAGACAGAGT  TTTTTTTTTTTTAGACAGAGT  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
    106985   . T   107006    |        0   106985  |  TTTTTTTTTT.TCCTGAGCAG  TTTTTTTTTTTTTCCTGAGCA  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
    106985   . T   107007    |        0   106985  |  TTTTTTTTTT.TCCTGAGCAG  TTTTTTTTTTTTCCTGAGCAG  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
    182920   . A   182943    |      564    47946  |  AAAAAAAAAA.ATGTGGTCTC  AAAAAAAAAAAATGTGGTCTC  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
    183484   . A   183508    |      564    47382  |  AAAAAAAAAA.ATAGATGAAC  AAAAAAAAAAAATAGATGAAC  |  1  1  CH17-157L1\tCH17-157L1|quiver (esc)
