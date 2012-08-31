

Variant Caller Functional Specification
=======================================

Version 1.4


Introduction
------------

This document describes the interface, input/output, and performance
characteristics of ``variantCaller.py``, a variant calling tool
provided by the ``GenomicConsensus`` package.


Software Overview
-----------------

The ``GenomicConsensus`` package provides a command-line tool,
``variantCaller.py``, which provides several variant-calling algorithms for
PacBio sequencing data.  ``variantCaller.py`` replaces ``EviCons`` and
``SmrtBayes``, the previous (haploid, diploid---respectively) variant callers
at PacBio.



Functional Requirements
-----------------------

Command-line interface
``````````````````````

``variantCaller.py`` is invoked from the command line.  For example, a simple
invocation is::

        variantCaller.py -j8 --algorithm=quiver \
                         -r lambdaNEB.fa        \
                         -o variants.gff        \
                         aligned_reads.cmp.h5   

which requests that variant calling proceed,
        - using 8 worker processes,
        - employing the **quiver** algorithm,
        - taking input from the file ``aligned_reads.cmp.h5``,
        - using the FASTA file ``lambdaNEB.fa`` as the reference,
        - and writing output to ``variants.gff``.

A particularly useful option is ``--referenceWindow/-w``: this option
allows the user to direct the tool to perform variant calling
exclusively on a *window* of the reference genome, where the


Invoking

:: 

    variantCaller.py --help

will provide a help message explaining all available options; they will be
documented here shortly.



Input and output
````````````````
``variantCaller.py`` requires two input files:

    - A file of reference-aligned reads in PacBio's standard cmp.h5 format;
    - A FASTA file that has been processed by ReferenceUploader.

The tool's output is formatted in the GFF format, as described in (how
to link to other file?).  External tools can be used to convert the
GFF file to a VCF or BED file---two other standard interchange formats
for variant calling.

.. note::

        **Input cmp.h5 file requirements**

        ``variantCaller.py`` requires its input cmp.h5 file to be
        be sorted.  An unsorted file can be sorting using the tool
        ``cmpH5Sort.py``.
        
        The *quiver* algorithm in ``variantCaller.py`` requires its
        input cmp.h5 file to have the following *pulse features*:
            - ``InsQV``,
            - ``SubsQV``,
            - ``DelQV``,
            - ``DelTag``,
            - ``MergeQV``.
        
        The *plurality* algorithm can be run on cmp.h5 files that lack
        these features.

The input file is the main argument to ``variantCaller.py``, while the output
file is provided as an argument to the ``-o`` flag.  For example,

::

        variantCaller.py aligned_reads.cmp.h5 -r lambda.fa  -o variants.gff

will read input from ``aligned_reads.cmp.h5``, using the reference
``lambda.fa``, and send output to the file ``variants.gff``.  The
extension of the filename provided to the ``-o`` flag is meaningful,
as it determines the output file format.  The file formats presently
supported, by extension, are

``.gff``
        GFFv3 format

``.txt``
        a simplified human readable format used primarily by the developers

If the ``-o`` flag is not provided, the default behavior is to output to a
``variants.gff`` in the current directory.


.. note::

    ``variantCaller.py`` does **not** modify its input cmp.h5 file
    in any way.  This is in contrast to previous variant callers in
    use at PacBio, which would write a *consensus* dataset to the input
    cmp.h5 file.


Available algorithms
````````````````````

At this time there are two algorithms available for variant calling:
**plurality** and **quiver**.  

**Plurality** is a simple and very fast procedure that merely tallies the most
frequent read base or bases found in alignment with each reference base, and
reports deviations from the reference as potential variants.

**Quiver** is a more complex procedure based on algorithms originally
developed for CCS.  Quiver leverages the quality values (QVs) provided by
upstream processing tools, which provide insight into whether
insertions/deletions/substitutions were deemed likely at a given read
position.  Use of **quiver** requires the ``ConsensusCore`` library as well as
trained parameter set, which will be loaded from a standard location (TBD).
Quiver can be thought of as a QV-aware local-realignment procedure.

Both algorithms are expected to converge to *zero* errors (miscalled variants)
as coverage increases; however **quiver** should converge much faster (i.e.,
fewer errors at low coverage), and should provide greater variant detection
power at a given error level.


Software interfaces
```````````````````
The ``GenomicConsensus`` module has two essential dependencies:

1. **pbcore**, the PacBio Python bioinformatics library
2. **ConsensusCore**, a C++ library with SWIG bindings that provides access to
   the same algorithms used in circular consensus sequencing.

Both of these modules are easily installed using their ``setup.py`` scripts,
which is the canonical means of installing Python packages.


Confidence values
-----------------

Both *quiver* and *plurality* make a confidence metric available for
every position of the consensus sequence.  The confidence should be
interpreted as a phred-transformed posterior probability that the
consensus call is incorrect; i.e.

.. math::

    QV = -10 \log_{10}(p_{err})

``variantCaller.py`` clips reported QV values at 93---larger values
cannot be encoded in a standard FASTQ file.


Performance Requirements
------------------------

``variantCaller.py`` performs variant calling in parallel using multiple
processes.  Work splitting and inter-process communication are handled using
the Python ``multiprocessing`` module.  Work can be split among an arbitrary
number of processes (using the ``-j`` command-line flag), but for best
performance one should use no more worker processes than there are CPUs in the
host computer.

The running time of the *plurality* algorithm should not exceed the
runtime of the BLASR process that produced the cmp.h5. The running
time of the *quiver* algorithm should not exceed 4x the runtime of
BLASR.

The amount of core memory (RAM) used among all the python processes launched
by a ``variantCaller.py`` run should not exceed the size of the uncompressed
input ``.cmp.h5`` file.




