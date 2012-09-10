
How-to install and use Quiver
=============================

Background
----------

*Quiver* is an algorithm for calling highly accurate consensus from
multiple PacBio raw reads, using a Hidden Markov Model exploiting
both the basecalls and QV metrics to infer the true underlying DNA
sequence.

*Quiver* is embodied within the `ConsensusCore` library, a C++ library
with C# and Python bindings.  To experiment with Quiver from Python,
you will need to install three packages in your Python environment:

- `pbcore`, a package providing access to PacBio data files
- `ConsensusCore`, containing Quiver
- `GenomicConsensus`, containing convenience routines and command-line tools

Step 1: Install required libraries
----------------------------------

The following are required:

- Boost  >= 1.4.7   (standard C++ libraries)
- SWIG   >= 2.0.7   (library wrapper generator)
- Python >= 2.7.2
- virtualenv        (builds isolated Python environments)

If you within PacBio, these requirements are already installed within
the cluster environment.

Otherwise, you will need to install them yourself.


Step 2: Set up your Python environment
--------------------------------------

I recommend using a Python *virtualenv* to isolate your sandbox.

To set up a new virtualenv, do ::

    $ cd; virtualenv -p python2.7 --no-site-packages VE-QUIVER

and activate the virtualenv using ::

    $ source ~/VE-QUIVER/bin/activate

There are some additional Python libraries required (NumPy and h5py),
which can be installed via ::

    $ pip install numpy==1.6.1
    $ pip install h5py==2.0.1

Step 3: Install PacBio libraries
--------------------------------

To install the PacBio software, execute ::

    $ pip install git+https://github.com/PacificBiosciences/pbcore
    $ git clone https://github.com/PacificBiosciences/ConsensusCore
    $ cd ConsensusCore; python setup.py install --swig=$SWIG --boost=$BOOST
    $ pip install git+https://github.com/PacificBiosciences/GenomicConsensus

where you replace $SWIG with the path to your SWIG executable and
$BOOST with the path to your boost install (the top level directory).
(Note that if SWIG is in your $PATH and boost is in `/usr/local` or
`/usr/include/`, you do not need to specify these flags on the command
line--`setup.py` will find them).


Step 4: Run Quiver
------------------
Experimental users are welcome to learn how to use the Quiver APIs by looking at

`GenomicConsensus/quiver/demo.py`

however note that this demo is optimized for didactic simplicity, not
consensus accuracy.  After understanding the demo code, look at
`quiver.py` to see how we handle edge cases that limit accuracy.

Those who wish to call consensus on a resequencing job can simply use
the `variantCaller.py` script that has been installed in your
virtualenv (from `GenomicConsensus`).

For example, ::

    $ variantCaller.py -j8 --algorithm=quiver         \
    >    aligned_reads.cmp.h5 -r path/to/lambda.fasta \
    >    -o variants.gff -o consensus.fasta

will use 8 CPUs to run Quiver on `aligned_reads.cmp.h5`, outputting
the consensus sequence and variants.

NOTE: the above command assumes your cmp.h5 was generated using the
Resequencing_QVs workflow in smrtportal.  If you have not run this
workflow, your cmp.h5 will not have the full complement of QV metrics
and therefore a different model will need to be used.  To use a model
that requires no QV metrics, try ::

    $ variantCaller.py -j8 --algorithm=quiver         \
    >    --parameters=NoQVsModel.trainedParams1       \
    >    aligned_reads.cmp.h5 -r path/to/lambda.fasta \
    >    -o variants.gff -o consensus.fasta

The QVs enable the highest accuracy consensus, especially in low
coverage regions, so we do not recommend using this



Step 5: Highly-accurate assembly consensus
------------------------------------------

Quiver enables consensus accuracies on genome assemblies at accuracies
greater than Q50 (one error per 100,000 bases).  At the present moment
this is not a streamlined workflow---we are working on this right now.

Here are the current steps needed to get a high-accuracy assembly
using Quiver:

- Start with a rough assembly---from any platform or assembly
  technique; PacBio users can use the `Allora` or `CeleraAssembler`
  workflows to generate a rough assembly.

- The output of the assembly is a FASTA file; at the present time the
  user must download this FASTA file and then import it as a new
  reference into SMRTPortal.

- Run a `Resequencing_QVs` job using the original data files, and the
  rough assembly FASTA file as a reference.

The output of the Resequencing_QVs job is the cmp.h5 file you will now
feed to Quiver::

    $ variantCaller.py -j8 --algorithm=quiver                 \
    >    aligned_reads.cmp.h5 -r path/to/rough-assembly.fasta \
    >    -o quiver-assembly.fasta

The `quiver-assembly.fasta` file contains the refined assembly. If you
have consisently high coverage across the genome, the quality should
be quite high.  Note that Quiver does *not* join contigs---it merely
refines their accuracy.

Step 5: Learn about Quiver
--------------------------

A presentation describing some of the details of how Quiver works is
available in `ConsensusCore/doc/Presentations/BrownBag2012/presentation.pdf`.
