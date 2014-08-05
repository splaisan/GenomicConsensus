
Run quiver on a large-insert C2 E. coli job.

  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/ecoli/job_059531.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/ecoli/ecoliK12_pbi_March2013.fasta
  $ quiver -j${JOBS-8} $INPUT -r $REFERENCE -o variants.gff -o css.fasta

Inspect the variants list.  A few mutations seem to have crept in
since I built the new reference.

  $ sed 's/\t/ /g' variants.gff
  ##gff-version 3
  ##pacbio-variant-version 2.1
  ##date * (glob)
  ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
  ##source GenomicConsensus * (glob)
  ##source-commandline * (glob)
  ##source-alignment-file * (glob)
  ##source-reference-file * (glob)
  ##sequence-region ecoliK12_pbi_March2013 1 4642522
  ecoliK12_pbi_March2013 . deletion 85 85 . . . reference=G;variantSeq=.;coverage=53;confidence=48
  ecoliK12_pbi_March2013 . deletion 219 219 . . . reference=A;variantSeq=.;coverage=58;confidence=47
  ecoliK12_pbi_March2013 . insertion 1536 1536 . . . reference=.;variantSeq=C;coverage=91;confidence=50

No no-call windows.

  $ fastacomposition css.fasta
  css.fasta A 1141540 C 1177642 G 1180362 T 1142977

MuMMer analysis.  No structural diffs

  $ nucmer -mum $REFERENCE css.fasta 2>/dev/null
  $ show-diff -H out.delta

SNPs same as variants

  $ show-snps -C -H out.delta  | sed  's/\s\+/ /g'
   85 G . 84 | 85 84 | 1 1 ecoliK12_pbi_March2013 ecoliK12_pbi_March2013|quiver
   220 A . 218 | 135 218 | 1 1 ecoliK12_pbi_March2013 ecoliK12_pbi_March2013|quiver
   1536 . C 1535 | 1316 1535 | 1 1 ecoliK12_pbi_March2013 ecoliK12_pbi_March2013|quiver
