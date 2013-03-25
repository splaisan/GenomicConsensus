
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
  ##sequence-region ecoliK12_mutated 1 4639560
  ecoliK12_mutated . substitution 547694 547694 . . . variantSeq=G;reference=A;coverage=100;confidence=48;length=1
  ecoliK12_mutated . insertion 547831 547831 . . . variantSeq=G;coverage=100;confidence=48;length=1
  ecoliK12_mutated . insertion 2171384 2171384 . . . variantSeq=CC;coverage=100;confidence=49;length=2
  ecoliK12_mutated . insertion 3422255 3422255 . . . variantSeq=C;coverage=65;confidence=51;length=1
  ecoliK12_mutated . insertion 4294288 4294288 . . . variantSeq=CG;coverage=100;confidence=48;length=2

OK, now let's take a gander at the consensus FASTA output.  We end up
with 17000 nocalls here, which is not bad considering the short insert
size (2KB) and the fact that E. Coli has 7 5KB near-perfect repeats.

  $ fastacomposition css.fasta
  css.fasta A 1137820 C 1175416 G 1172314 T 1137029 a 4381 c 4115 g 4586 t 3918

Use the MuMMer suite to look at the differences from the reference.

  $ nucmer -mum $REFERENCE css.fasta 2>/dev/null

First: no structural differences.  There are some (gaps) if you use
the masked output.

  $ show-diff -q out.delta | sed 's/\t/ /g'
  * (glob)
  NUCMER
  
  [SEQ] [TYPE] [S1] [E1] [LEN 1]

Next, the SNPs.

  $ show-snps -H -C -x10 out.delta
    547694   A G   547692    |      140   547692  |  AACTAAACGAAGGCAAAACAC  AACTAAACGAGGGCAAAACAC  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
    547834   . G   547833    |      140   547833  |  TGGAGCAGGG.GAAGTGAACT  TGGAGCAGGGGGAAGTGAACT  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   2171385   . C   2171385   |        0  2171385  |  CCGATACCAC.CGCCGTATGT  CCGATACCACCCCGCCGTATG  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   2171385   . C   2171386   |        0  2171385  |  CCGATACCAC.CGCCGTATGT  CGATACCACCCCGCCGTATGT  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   3619034   . C   3619040   |   193936  1020527  |  GTCATTGCCC.CGGACGGCAG  GTCATTGCCCCCGGACGGCAG  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   4209525   . C   4209540   |       13   430036  |  ACGGTTGTCC.CGGTTTAAGC  ACGGTTGTCCCCGGTTTAAGC  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   4209538   . T   4209554   |       13   430023  |  TTTAAGCGTG.TAGGCTGGTT  TTTAAGCGTGTTAGGCTGGTT  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   4209580   . C   4209597   |       42   429981  |  AAGGCTGAGG.CGTGATGACG  AAGGCTGAGGCCGTGATGACG  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   4294288   . C   4294306   |        0   345273  |  AAGGCGTTTA.CCGCATCCGA  AAGGCGTTTACGCCGCATCCG  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
   4294288   . G   4294307   |        0   345273  |  AAGGCGTTTA.CCGCATCCGA  AGGCGTTTACGCCGCATCCGA  |  1  1  ecoliK12_mutated\tecoliK12_mutated|quiver (esc)
