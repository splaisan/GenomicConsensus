
  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/ecoli/job_044601.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/ecoli/ecoli_mutated.fasta
  $ quiver -j${JOBS-8} $INPUT -r $REFERENCE -o variants.gff -o css.fasta

NOTE: We use sed to replace tabs with spaces before inspecting output
in some of the tests below.  It just makes the test file more legible.

Inspect the variants list.  At the moment we are left with a short
list of Quiver errors in dinucleotide repeat regions (fix will be
brought into mainline shortly), and true variants (errors in the
reference!).

  $ sed 's/\t/ /g' variants.gff
  ##gff-version 3
  ##pacbio-variant-version 1.4
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus * (glob)
  ##source-commandline * (glob)
  ##sequence-region ref000001|ecoliK12_mutated 1 4639560
  ref000001|ecoliK12_mutated . substitution 547694 547694 . . . variantSeq=G;reference=A;coverage=100;confidence=48;length=1
  ref000001|ecoliK12_mutated . insertion 547831 547831 . . . variantSeq=G;coverage=100;confidence=48;length=1
  ref000001|ecoliK12_mutated . insertion 2171384 2171384 . . . variantSeq=CC;coverage=100;confidence=49;length=2
  ref000001|ecoliK12_mutated . insertion 3422255 3422255 . . . variantSeq=C;coverage=65;confidence=51;length=1
  ref000001|ecoliK12_mutated . insertion 4294288 4294288 . . . variantSeq=CG;coverage=100;confidence=48;length=2

OK, now let's take a gander at the consensus FASTA output.  We end up
with 17000 nocalls here, which is not bad considering the short insert
size (2KB) and the fact that E. Coli has 7 5KB near-perfect repeats.

  $ fastacomposition css.fasta
  css.fasta A 1137820 C 1175416 G 1172314 N 17000 T 1137029

Use the MuMMer suite to look at the differences from the reference.

  $ nucmer -mum $REFERENCE css.fasta 2>/dev/null

First, structural differences.  Not really sure why some of the Ns are
deemed JMP vs GAP here.

  $ show-diff -q out.delta | sed 's/\t/ /g'
  * (glob)
  NUCMER

  [SEQ] [TYPE] [S1] [E1] [LEN 1]
  ref000001|ecoliK12_mutated|quiver JMP 225001 226000 1000
  ref000001|ecoliK12_mutated|quiver GAP 226500 226999 500 500 0
  ref000001|ecoliK12_mutated|quiver GAP 257499 257998 500 500 0
  ref000001|ecoliK12_mutated|quiver GAP 1298500 1298999 500 500 0
  ref000001|ecoliK12_mutated|quiver GAP 1871000 1871499 500 500 0
  ref000001|ecoliK12_mutated|quiver JMP 2725502 2727501 2000
  ref000001|ecoliK12_mutated|quiver JMP 3423007 3425006 2000
  ref000001|ecoliK12_mutated|quiver GAP 3618506 3619005 500 500 0
  ref000001|ecoliK12_mutated|quiver GAP 3761507 3762506 1000 1000 0
  ref000001|ecoliK12_mutated|quiver JMP 3941007 3943006 2000
  ref000001|ecoliK12_mutated|quiver GAP 4034509 4036508 2000 2000 0
  ref000001|ecoliK12_mutated|quiver GAP 4037009 4037508 500 500 0
  ref000001|ecoliK12_mutated|quiver JMP 4166010 4168009 2000
  ref000001|ecoliK12_mutated|quiver JMP 4207515 4209514 2000


Next, the SNPs.  Most of these are low confidence SNPs due to coverage
dropouts.  Just count them for now.

  $ show-snps -H -C -x10 out.delta | wc -l
  20
