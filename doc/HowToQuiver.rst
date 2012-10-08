
How to install and use Quiver
=============================

*Note: please install this software on an isolated machine that does
not have SMRTanalysis installed.  Older versions of SMRTanalysis
pollute the ``PYTHONPATH``, which has the undesirable effect of
overriding ``virtualenv``-installed modules.*

Background
----------
**Quiver** is an algorithm for calling highly accurate consensus from
multiple PacBio raw reads, using a Hidden Markov Model exploiting both
the basecalls and QV metrics to infer the true underlying DNA
sequence.

Quiver is available through the ``variantCaller.py`` script from the
``GenomicConsensus`` package, using the ``--algorithm=quiver`` option.
To use Quiver, the following PacBio software is required.

- ``GenomicConsensus``, containing ``variantCaller.py``
- ``ConsensusCore``, a C++ library containing the core computational
  routines for Quiver
- ``pbcore``, a package providing access to PacBio data files


Required libraries and tools
----------------------------
To install the PacBio software, the following are required:

- Boost  >= 1.4.7   (standard C++ libraries)
- SWIG   >= 2.0.7   (library wrapper generator)
- Python >= 2.7.2
- virtualenv        (builds isolated Python environments)

If you within PacBio, these requirements are already installed within
the cluster environment.

Otherwise, you will need to install them yourself.  The automatic
installation script requires that the ``swig`` executable is in your
UNIX ``$PATH`` and that your boost installation can be found under
``/usr/include`` or ``/usr/local``.


Data file requirements
----------------------

To make the most accurate consensus calls possible, Quiver makes use
of a battery of quality value metrics calculated by the basecaller.
By default, only a subset of these quality values are included in the
``.cmp.h5`` files produced by SMRTanalysis.  To get a ``.cmp.h5`` with
all the QVs loaded, you will need to use the ``Resequencing_QVs``
workflow within SMRTPortal.

If you are using SMRTPortal 1.3.3 or later, ``Resequencing_QVs`` will
already be installed.

Otherwise, it is easy to install into your SMRTPortal installation.
All you need to do is put the `Resequencing_QVs.1.xml`_ file into
SMRTportal's ``common/protocols/`` directory, and `BLASR_QVs.1.xml`_ in
the ``common/protocols/mapping`` directory.

Automatic installation instructions
-----------------------------------
If your system meets the installation requirements, you can perform an
automatic installation of the PacBio software for Quiver by
executing::

    $ curl -L http://git.io/JR7TnQ | bash


Manual installation instructions
--------------------------------
If your SWIG or BOOST installations are in non-standard locations or
you encounter a problem with the automatic installation script, you
can follow these steps to install Quiver manually.



Step 1: Set up your Python environment
``````````````````````````````````````
I recommend using a Python *virtualenv* to isolate your sandbox.

To set up a new virtualenv, do ::

    $ cd; virtualenv -p python2.7 --no-site-packages VE-QUIVER

and activate the virtualenv using ::

    $ source ~/VE-QUIVER/bin/activate

There are some additional Python libraries required (NumPy and h5py),
which can be installed via ::

    $ pip install numpy==1.6.1
    $ pip install h5py==2.0.1


Step 2: Install PacBio libraries
````````````````````````````````
To install the PacBio software, execute ::

    $ pip install git+https://github.com/PacificBiosciences/pbcore
    $ git clone https://github.com/PacificBiosciences/ConsensusCore
    $ cd ConsensusCore; python setup.py install --swig=$SWIG --boost=$BOOST
    $ pip install git+https://github.com/PacificBiosciences/GenomicConsensus

where you replace ``$SWIG`` with the path to your ``swig`` executable
and ``$BOOST`` with the path to your boost install (the top level
directory).  (Note that if SWIG is in your ``$PATH`` and boost is in
``/usr/local`` or ``/usr/include/``, you do not need to specify these
flags on the command line---``setup.py`` will find them).


Step 3: Run Quiver
``````````````````
Those who wish to call consensus on a resequencing job can simply use
the ``variantCaller.py`` script that has been installed in your
virtualenv (from `GenomicConsensus`).

For example, ::

    $ variantCaller.py -j8 --algorithm=quiver         \
    >    aligned_reads.cmp.h5 -r path/to/lambda.fasta \
    >    -o variants.gff -o consensus.fasta

will use 8 CPUs to run Quiver on ``aligned_reads.cmp.h5``, outputting
the consensus sequence and variants.

Note that if you have not used the `Resequencing_QVs` protocol to
generate the cmp.h5 file---or if the source bas.h5 file was generated
by pre-1.3.1 instrument software---the cmp.h5 will not contain the
full battery of QV metrics required for optimal Quiver accuracy.  the
command will still work, but it will give a warning that its accuracy
will be suboptimal.


Step 4: Highly-accurate assembly consensus
``````````````````````````````````````````
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

- Run a `Resequencing_QVs` job using the original data files, and
  the rough assembly FASTA file as a reference.

The output of the `Resequencing_QVs` job is the cmp.h5 file you will now
feed to Quiver::

    $ variantCaller.py -j8 --algorithm=quiver                 \
    >    aligned_reads.cmp.h5 -r path/to/rough-assembly.fasta \
    >    -o quiver-assembly.fasta

The ``quiver-assembly.fasta`` file contains the refined assembly. If
you have consisently high coverage across the genome, the quality
should be quite high.  Note that Quiver does *not* join contigs---it
merely refines their accuracy.


Learn About Quiver
------------------

There is currently a presentation_ available describing some of the
details of how Quiver works.  More documentation is in the works.

Experimental users are welcome to learn how to use the Quiver APIs by
read the source file ``GenomicConsensus/quiver/demo.py``

However, note that this demo is optimized for didactic simplicity, not
consensus accuracy.  After understanding the demo code, look at
``quiver.py`` to see how we handle edge cases that limit accuracy.


.. _presentation: https://github.com/PacificBiosciences/ConsensusCore/raw/master/doc/Presentations/BrownBag2012/presentation.pdf
.. _Resequencing_QVs.1.xml: https://github.com/PacificBiosciences/ConsensusCore/raw/master/tools/Resequencing_QVs.1.xml
.. _BLASR_QVs.1.xml: https://github.com/PacificBiosciences/ConsensusCore/raw/master/tools/BLASR_QVs.1.xml


Known Issues
------------

``ConsensusCore`` doesn't compile with Boost 1.48.  1.47 and 1.51 are
fine.
