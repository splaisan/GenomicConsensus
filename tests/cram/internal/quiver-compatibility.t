
Quiver should abort if the cmp.h5 is not suitable.  Let's make sure it does the right thing.

First make sure we abort once we recognize the tiny fluidigm file is CCS data.

  $ export DATA=$TESTDIR/../../data
  $ export INPUT=$DATA/fluidigm_amplicons/040500.cmp.h5
  $ export REFERENCE=$DATA/fluidigm_amplicons/Fluidigm_human_amplicons.fasta
  $ quiver -r $REFERENCE -o variants.gff $INPUT 2>1
  Failure: The Quiver algorithm requires a cmp.h5 file containing standard (non-CCS) reads.
  [255]


Tiny lambda file.  Make sure it recognizes this cmp.h5 has an imcomplete set of QVs.

  $ export INPUT=/mnt/secondary/Share/Quiver/TestData/tinyLambda/aligned_reads_1.cmp.h5
  $ export REFERENCE=/mnt/secondary/Share/Quiver/TestData/tinyLambda/lambdaNEB.fa
  $ quiver -p C2.AllQVsModel -r $REFERENCE -o variants.gff $INPUT
  Failure: Selected Quiver parameter set is incompatible with this cmp.h5 file due to missing data tracks.
  [255]


It should handle the request of a parameter set by complete name:

  $ quiver -v -p C2.NoQVsModel -r $REFERENCE -o variants.gff $INPUT 2>&1 | grep "Using Quiver parameter set"
  [INFO] Using Quiver parameter set(s) C2.NoQVsModel

... or by chemistry name:

  $ quiver -v -p C2 -r $REFERENCE -o variants.gff $INPUT 2>&1 | grep "Using Quiver parameter set"
  [INFO] Using Quiver parameter set(s) C2.NoQVsModel


... and should fail informatively when we ask for an unrecognized
parameter set or chemistry:

  $ quiver -p SuperChem.Model -r $REFERENCE -o variants.gff $INPUT
  Quiver: no available parameter set named SuperChem.Model
  [255]

  $ quiver -p SuperChem -r $REFERENCE -o variants.gff $INPUT
  Quiver: no parameter set available compatible with this cmp.h5 for chemistry "SuperChem" 
  [255]
