

Variant Caller Functional Specification
=======================================

Version 3.3


Introduction
------------

This document describes the interface, input/output, and performance
characteristics of ``variantCaller``, a variant calling tool
provided by the ``GenomicConsensus`` package.


Software Overview
-----------------

The ``GenomicConsensus`` package provides a command-line tool,
``variantCaller``, which provides several consensus and variant-calling  algorithms for
PacBio sequencing data.


Functional Requirements
-----------------------

Command-line interface
``````````````````````

``variantCaller`` is invoked from the command line.  For example, a simple
invocation is::

        variantCaller -j8 --algorithm=arrow  \
                         -r lambdaNEB.fa     \
                         -o variants.gff     \
                         aligned_subreads.bam

which requests that variant calling proceed,
- using 8 worker processes,
- employing the **arrow** algorithm,
- taking input from the file ``aligned_subreads.bam``,
- using the FASTA file ``lambdaNEB.fa`` as the reference,
- and writing output to ``variants.gff``.

A particularly useful option is ``--referenceWindow/-w``: this option
allows the user to direct the tool to perform variant calling
exclusively on a *window* of the reference genome.

Invoking

::

    variantCaller --help

will provide a help message explaining all available options; they will be
documented here shortly.



Input and output
````````````````
``variantCaller`` requires two input files:

- A sorted file of reference-aligned reads in `PacBio's standard BAM format`_;
- A FASTA file adhering to `PacBio's FASTA file conventions`_

The input file is the main argument to ``variantCaller``, while the
output files are provided as arguments to the ``-o`` flag.  For
example,

::

        variantCaller aligned_subreads.bam -r lambda.fa  -o myVariants.gff -o myConsensus.fasta

will read input from ``aligned_subreads.bam``, using the reference
``lambda.fa``, and send variant call output to the file
``myVariants.gff``, and consensus output to ``myConsensus.fasta``.
The extension of the filename provided to the ``-o`` flag is
meaningful, as it determines the output file format.  The file formats
presently supported, by extension, are

``.gff``
        PacBio GFFv3 variants format; convertable to VCF or BED.

``.fasta``
        FASTA file recording the consensus sequence calculated for each reference contig

``.fastq``
        FASTQ file recording the consensus sequence calculated for
        each reference contig, as well as per-base confidence scores


.. note::

   The *quiver* and *arrow* algorithms require that certain metrics
   are in place in the input BAM file.

   *quiver*, which operates on RSII data only, requires the
   basecaller-computed "pulse features" ``InsertionQV``,
   ``SubstitutionQV``, ``DeletionQV``, and ``DeletionTag``.  These
   features are populated in BAM tags by the ``bax2bam`` conversion
   program.

   *arrow*, which operates on RSII P6-C4 data and all Sequel data,
   requires per-read SNR metrics, and the per-base ``PulseWidth``
   metric for Sequel data (but not for RSII P6-C4).  These metrics are
   populated by Sequel instrument software or the ``bax2bam``
   converter (for RSII data).

   The selected algorithm will halt with an error message if features
   it requires are unavailable.


Available algorithms
````````````````````

At this time there are two algorithms available for variant calling:
**plurality**, **quiver**, and **arrow**.

**Plurality** is a simple and very fast procedure that merely tallies
the most frequent read base or bases found in alignment with each
reference base, and reports deviations from the reference as potential
variants.  This is a very insensitive and flawed approach for PacBio
sequence data, which is prone to insertion and deletion errors.

**Quiver** is a more complex procedure based on algorithms originally
developed for CCS.  Quiver leverages the quality values (QVs) provided by
upstream processing tools, which provide insight into whether
insertions/deletions/substitutions were deemed likely at a given read
position.  Use of **quiver** requires the ``ConsensusCore``
library.

**Arrow** is the successor to Quiver; it uses a more principled HMM
model approach.  It does not require basecaller quality value metrics;
rather, it uses the per-read SNR metric and the per-pulse pulsewidth
metric as part of its likelihood model.  Beyond the model specifics,
other aspects of the Arrow algorithm are similar to Quiver.  Use of
**arrow** requires the ``ConsensusCore2`` library, which is provided
by the ``unanimity`` codebase.


Software interfaces
```````````````````
The ``GenomicConsensus`` module has two essential dependencies:

1. **pbcore**, the PacBio Python bioinformatics library
2. **ConsensusCore**, a C++ library with SWIG bindings that provides
   access to the Quiver algorithm.
3. **ConsensusCore2**, a C++ library with SWIG bindings that provides access to
   the Arrow algorithm.

These modules are easily installed using their ``setup`` scripts,
which is the canonical means of installing Python packages.


Confidence values
-----------------

The arrow*, *quiver*, and *plurality* algorithms make a confidence
metric available for every position of the consensus sequence.  The
confidence should be interpreted as a phred-transformed posterior
probability that the consensus call is incorrect; i.e.

.. math::

    QV = -10 \log_{10}(p_{err})

``variantCaller`` clips reported QV values at 93---larger values
cannot be encoded in a standard FASTQ file.



Chemistry specificity
---------------------

The Quiver and Arrow algorithm parameters are trained per-chemistry.
Quiver and Arrow identify the sequencing chemistry used for each run
by looking at metadata contained in the data file (the input BAM or
cmp.h5 file).  This behavior can be overriden by a command line flag.

When multiple chemistries are represented in the reads in the input
file, Quiver/Arrow will model each read appropriately using the
parameter set for its chemistry, thus yielding optimal results.


Performance Requirements
------------------------

``variantCaller`` performs variant calling in parallel using multiple
processes.  Work splitting and inter-process communication are handled using
the Python ``multiprocessing`` module.  Work can be split among an arbitrary
number of processes (using the ``-j`` command-line flag), but for best
performance one should use no more worker processes than there are CPUs in the
host computer.

The running time of the *plurality* algorithm should not exceed the
runtime of the BLASR process that produced the cmp.h5. The running
time of the *quiver* algorithm should not exceed 4x the runtime of
BLASR.

The amount of core memory (RAM) used by a ``variantCaller`` run should
not exceed 2GB per active CPU core (as selected using the ``-j`` flag).


.. _PacBio's standard BAM format: http://pacbiofileformats.readthedocs.io/en/3.0/BAM.html
.. _PacBio's FASTA file conventions: http://pacbiofileformats.readthedocs.io/en/3.0/FASTA.html
