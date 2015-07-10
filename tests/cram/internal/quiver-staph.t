This input data is taken from the output of the mapping job in Pysiv:
pysiv_jobs/jobs/BAMMapping/saureus_p6c4

  $ export BAM=/mnt/secondary/Share/Quiver/TestData/staph/m140911_084715_42139_c100702390480000001823141103261514_s1_p0.aligned_subreads.bam
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/staph/S_aureus_USA300_TCH1516.fasta
  $ export MASK=/mnt/secondary/Share/Quiver/GenomeMasks/S_aureus_USA300_TCH1516-mask.gff 

  $ quiver -j${JOBS-8} $BAM -r $REFERENCE -o variants.gff -o css.fasta -o css.fastq

Inspect the variant calls.  The first variant call might be an error
(follow up on this) but the latter is an error in the reference, it
seems.

  $ gffsubtract.pl variants.gff $MASK | grep -v "#" | sed 's/\t/ /g'

One window is nocalled.  Follow up on this.

  $ fastacomposition css.fasta
  css.fasta A 960306 C 470724 G 470270 T 971459

One gap

  $ nucmer -mum $REFERENCE css.fasta 2>/dev/null
  $ show-diff -H out.delta | sed 's/\t/ /g'
  Staphylococcus_aureus_subsp_aureus_USA300_TCH1516 GAP 2149110 2149328 219 67 152

(Not showing the SNPs here as they just correspond to the masked-out region that we know doesn't match the reference)
