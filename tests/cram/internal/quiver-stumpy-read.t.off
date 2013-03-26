
We occasionally come across a "stumpy" read, where there is a large
gap in the read, due to a single molecule event that is as-yet not
well understood.  At the moment if these guys aren't identified and
filtered out the POA and probably the rest of Quiver will crash and
burn.  This is a simple test to make sure we get it right.  The file
contains coverage restricted to a small ~10KB window containing the
stumpy read.

  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/stumpyReadInEcoli/out.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/stumpyReadInEcoli/ecoli_mutated.fasta
  $ quiver --noEvidenceConsensusCall=nocall -j${JOBS-8} $INPUT -r $REFERENCE \
  > -o variants.gff -o css.fasta

Now compare back to the reference.  Coverage island results in
consensus sequence only in the window [2734568,2745424] (GFF
convention).

  $ nucmer -mum $REFERENCE css.fasta 2>/dev/null
  $ show-coords -H out.delta
   2734568  2745424  |  2734568  2745424  |    10857    10857  |    99.96  | ecoliK12_mutated\tecoliK12_mutated|quiver (esc)

No confident variants.  The 4 SNPs in the alignment are at the
fringes, where there is low coverage.

  $ grep -v "#" variants.gff
  [1]
