
This input data is taken from the output of the mapping job in Pysiv:
pysiv_jobs/jobs/BAMMapping/saureus_p6c4

  $ export BAM=/mnt/secondary/Share/Quiver/TestData/staph/aligned_reads.bam
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/staph/S_aureus_USA300_TCH1516.fasta

  $ quiver -j${JOBS-8} $BAM -r $REFERENCE -o variants.gff -o css.fasta -o css.fastq
  [WARNING] 'fancyChunking' not yet available for BAM, disabling

Inspect the variant calls.  The first variant call might be an error
(follow up on this) but the latter is an error in the reference, it
seems.

  $ grep -v "#" variants.gff | sed 's/\t/ /g'
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 1592553 1592553 . . . reference=T;variantSeq=.;coverage=100;confidence=50
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 . deletion 2179435 2179435 . . . reference=C;variantSeq=.;coverage=100;confidence=49

One window is nocalled.  Follow up on this.

  $ fastacomposition css.fasta
  css.fasta A 960147 C 470678 G 470218 T 971370 a 176 c 74 g 91 t 159

No gaps.

  $ nucmer -mum $REFERENCE css.fasta 2>/dev/null
  $ show-diff -H out.delta | sed 's/\t/ /g'

Why is the latter variant not included here, below?  Is it because the
alignment is ambiguous?

  $ show-snps -H -C out.delta
   1592556   T .   1592555   |   586882  1280359  |  1  1  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516\tStaphylococcus_aureus_subsp_aureus_USA300_TCH1516|quiver (esc)
