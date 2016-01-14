
Test the (augmentation) of the alignment_summary.gff file by summarizeConsensus

  $ export DATA=/mnt/secondary/Share/Quiver/TestData/tinyLambda/
  $ export PATH=$TESTDIR/..:$PATH
  $ export VARIANTSGFF=$DATA/variants.gff.gz
  $ export ALIGNMENTSUMMARYGFF=$DATA/alignment_summary.gff
  $ summarizeConsensus              \
  >   --variantsGff $VARIANTSGFF    \
  >   $ALIGNMENTSUMMARYGFF          \
  >   -o alignment_summary.out.gff

  $ head -20 alignment_summary.out.gff
  ##gff-version 3
  ##date Thu, 03-Feb-2011 14:54:12
  ##source PACBIO_AlignmentSummary 1.0
  ##source ConsensusStats v0.1
  ##source-commandline summarizeCoverage.py --reference /mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/references/lambda --numRegions=500 /mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/jobs/016/016789/data/aligned_reads.cmp.h5
  ##source-commandline mono ConsensusStats.exe /mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/jobs/016/016789/data/variants.gff /mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/jobs/016/016789/data/aligned_reads.cmp.h5
  ##sequence-region lambda_NEB3011 1 48502
  ##source GenomicConsensus * (glob)
  ##pacbio-alignment-summary-version 0.6
  ##source-commandline * (glob)
  lambda_NEB3011\t.\tregion\t1\t100\t0.00\t+\t.\tcov2=150.440,26.772;gaps=0,0;cov=51,160,171;cQv=20,20,20;del=0;ins=19;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t101\t200\t0.00\t+\t.\tcov2=168.700,1.780;gaps=0,0;cov=166,168,173;cQv=20,20,20;del=0;ins=16;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t201\t300\t0.00\t+\t.\tcov2=167.860,1.732;gaps=0,0;cov=165,168,171;cQv=20,20,20;del=1;ins=17;sub=1 (esc)
  lambda_NEB3011\t.\tregion\t301\t400\t0.00\t+\t.\tcov2=177.690,2.587;gaps=0,0;cov=168,179,181;cQv=20,20,20;del=2;ins=8;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t401\t500\t0.00\t+\t.\tcov2=179.730,1.248;gaps=0,0;cov=177,180,182;cQv=20,20,20;del=0;ins=0;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t501\t600\t0.00\t+\t.\tcov2=186.670,4.907;gaps=0,0;cov=177,188,195;cQv=20,20,20;del=0;ins=0;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t601\t700\t0.00\t+\t.\tcov2=200.160,4.051;gaps=0,0;cov=192,200,206;cQv=20,20,20;del=0;ins=0;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t701\t800\t0.00\t+\t.\tcov2=213.630,7.634;gaps=0,0;cov=200,215,226;cQv=20,20,20;del=0;ins=0;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t801\t900\t0.00\t+\t.\tcov2=244.290,12.954;gaps=0,0;cov=224,243,262;cQv=20,20,20;del=0;ins=0;sub=0 (esc)
  lambda_NEB3011\t.\tregion\t901\t1000\t0.00\t+\t.\tcov2=267.070,3.724;gaps=0,0;cov=259,266,274;cQv=20,20,20;del=0;ins=0;sub=0 (esc)

  $ grep -v '\#.*' alignment_summary.out.gff | md5sum
  08f89b262b159671cdcdd8bdc8331461  -
