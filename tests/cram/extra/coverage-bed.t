Test conversion of alignment summary GFF to coverage BED.

  $ export DATA=$TESTDIR/../../data
  $ export INPUT=$DATA/fluidigm_amplicons/alignment_summary.gff

  $ gffToBed --name=coverage \
  >          --description="PacBio coverage" \
  > coverage $INPUT > coverage.bed
  $ head -20 coverage.bed
  track name=coverage description="PacBio coverage" useScore=0
  ref000001\t0\t1\tmeanCov\t27.000\t+ (esc)
  ref000001\t1\t2\tmeanCov\t27.000\t+ (esc)
  ref000001\t2\t3\tmeanCov\t27.000\t+ (esc)
  ref000001\t3\t4\tmeanCov\t27.000\t+ (esc)
  ref000001\t4\t5\tmeanCov\t27.000\t+ (esc)
  ref000001\t5\t6\tmeanCov\t27.000\t+ (esc)
  ref000001\t6\t7\tmeanCov\t27.000\t+ (esc)
  ref000001\t7\t8\tmeanCov\t27.000\t+ (esc)
  ref000001\t8\t9\tmeanCov\t27.000\t+ (esc)
  ref000001\t9\t10\tmeanCov\t27.000\t+ (esc)
  ref000001\t10\t11\tmeanCov\t27.000\t+ (esc)
  ref000001\t11\t12\tmeanCov\t27.000\t+ (esc)
  ref000001\t12\t13\tmeanCov\t27.000\t+ (esc)
  ref000001\t13\t14\tmeanCov\t27.000\t+ (esc)
  ref000001\t14\t15\tmeanCov\t27.000\t+ (esc)
  ref000001\t15\t16\tmeanCov\t27.000\t+ (esc)
  ref000001\t16\t17\tmeanCov\t27.000\t+ (esc)
  ref000001\t17\t18\tmeanCov\t27.000\t+ (esc)
  ref000001\t18\t19\tmeanCov\t27.000\t+ (esc)
