
Quiver should abort if the cmp.h5 is not suitable.  Let's make sure it does the right thing.

First, make sure it recognizes this cmp.h5 has an imcomplete set of QVs.

  $ export DATA=$TESTDIR/../data
  $ export INPUT=$DATA/lambda/aligned_reads_1.cmp.h5
  $ export REFERENCE=$DATA/lambda/lambdaNEB.fa
  $ variantCaller.py --parameters=AllQVsModel.C2 --algorithm=quiver -r $REFERENCE -o variants.gff $INPUT
  Failure: CmpH5 file is incompatible with algorithm "Quiver": This Quiver parameter set requires QV features not available in this .cmp.h5 file.
  [255]

Next, make sure we abort once we recognize this file is CCS data.

  $ export INPUT=$DATA/fluidigm_amplicons/040500.cmp.h5
  $ export REFERENCE=$DATA/fluidigm_amplicons/Fluidigm_human_amplicons.fasta
  $ variantCaller.py --algorithm=quiver -r $REFERENCE -o variants.gff $INPUT 2>1
  Failure: CmpH5 file is incompatible with algorithm "Quiver": The Quiver algorithm requires a cmp.h5 file containing standard (non-CCS) reads.
  [255]
