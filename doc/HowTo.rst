
How to install and use GenomicConsensus
=======================================

**We recommend that you obtain GenomicConsensus by installing the most
recent version of SMRTanalysis.  Other means of installation are not
officially supported.**


Basic running instructions
--------------------------

Basic usage---using 8 CPUs to compute consensus of mapped reads and
variants relative to a reference---is as follows::

    % quiver -j8 aligned_reads{.cmp.h5, .bam, .fofn, or .xml} \
    >     -r reference{.fasta or .xml} -o variants.gff        \
    >     -o consensus.fasta -o consensus.fastq

``quiver`` is a shortcut for ``variantCaller --algorithm=quiver``.
Naturally, to use arrow you could use the ``arrow`` shortcut or
``variantCaller --algorithm=arrow``.

in this example we perform haploid consensus and variant calling on
the mapped reads in the ``aligned_reads.bam`` which was aligned to
``reference.fasta``.  The ``reference.fasta`` is only used for
designating variant calls, not for computing the consensus.  The
consensus quality score for every position can be found in the output
FASTQ file.

*Note that 2.3 SMRTanalysis does not support "dataset" input (FOFN
or XML files); those who need this feature should wait for the forthcoming
release of SMRTanalysis 3.0 or build from GitHub sources.*


Running a large-scale resequencing/polishing job in SMRTanalysis 2.3
--------------------------------------------------------------------

We do not recommend attempting  to construct a single giant cmp.h5 file and
then processing it on a single node.  This is inefficient and users attempting to do this
have run into many problems with the instability of the HDF5 library (which PacBio is
moving away from, in favor of BAM_.)

To run a large-scale resequencing job (>50 megabase genome @ 50x
coverage,nominally), you want to spread the computation load across
multiple nodes in your computing cluster.  

The `smrtpipe` workflow engine in SMRTanalysis 2.3 provides a
convenient workflow automating this---it will automatically spread the
load for both mapping and quiver jobs among your available cluster
nodes.  This is accessible via the SMRTportal UI; the simplest way to 
set up and run thse workflows is via tha UI.  Nonetheless, we include 
command-line instructions for completeness.

If you have to run the `smrtpipe` workflow manually from the command
line, a recipe is as folows::

1. Make sure the reference you will align and compare against is
   present in a SMRTportal "reference repository".  Even if you
   don't want to use SMRTportal, you need to build/import the
   reference appropriately, and the simplest way to do that is
   via SMRTportal.  If you don't have a SMRTportal instance, 
   you can use the ``referenceUploader`` command to prepare your
   reference repository.

2. Prepare an "input.fofn" file listing, one-per-line, each "bax.h5"
   file in your input data set.

3. Convert the "input.fofn" to an "input.xml" file that SMRTpipe can
   understand::

   $ fofnToSmrtpipeInput.py input.fofn > input.xml

4. Prepare your "params.xml" file.  Here is a `params.xml template`_
   you can use; you should just need to edit the reference path.

5. Activate your SMRTanalysis environment, and invoke smrtpipe::

   $ source <SMRT Analysis>/etc/setup.sh
   $ smrtpipe.py --distribute --params=params.xml xml:input.xml

6. After successful execution is complete, the results should be
   available as `data/consensus.fast[aq].gz` and
   `data/variants.gff.gz`, etc.

Please consult the `SMRTpipe reference manual`_ for further information.

*Note that resequencing (mapping reads against a reference genome and
then calling consensus and identifying variants) and polishing
(mapping reads against a draft assembly and then taking the consensus
output as the final, polished, assembly) are the same algorithmic
operation, the only effective difference is that the "variants.gff"
output is not biologically meaningful in the polishing case---it just
records the edits that were made to the draft to produce the polished
assembly.*

Running a large-scale quiver/arrow job in SMRTanalysis 3.0+
-----------------------------------------------------------

(Forthcoming)


Building bleeding-edge code (unsupported)
----------------------------------------

If you need to access the the latest code for some reason, a
convenient way to build it is to use PacBio's pitchfork_ build
system, which will take care of all third party dependencies for you.
Here's a recipe::

  git clone git@github.com:PacificBiosciences/pitchfork.git
  cd pitchfork
  make GenomicConsensus   # may take some time, as it builds dependencies...

Now, with GenomicConsensus built, you can use it via::

  bash --init-file deployment/setup-env.sh  # Puts you in a subshell where your build is available
  quiver --help                             # now you have quiver, arrow, etc. available

If you encounter build issues using `pitchfork`, please report the
issues there.  Note that you can deploy PacBio software to a location
of your choice using pitchfork.


Further questions?
------------------

Please consult the `FAQ document`_.

.. _`FAQ document`: ./FAQ.rst
.. _pitchfork : https://github.com/PacificBiosciences/pitchfork
.. _`params.xml template`: ./params-template.xml
.. _`SMRTpipe reference manual`: http://www.pacb.com/wp-content/uploads/2015/09/SMRT-Pipe-Reference-Guide.pdf
.. _`BAM`: http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
