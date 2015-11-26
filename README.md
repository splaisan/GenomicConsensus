GenomicConsensus (quiver)
-------------------------

The ``GenomicConsensus`` package provides the ``quiver`` tool,
PacBio's flagship consensus and variant caller.  The backend logic is
provided by the ``ConsensusCore`` and ``ConsensusCore2`` libraries,
which you must install first.


Installing
----------
Make sure you have set up and activated your virtualenv, and installed
``pbcore``, ``ConsensusCore``, and ``ConsensusCore2`` (which cannot be
installed automatically by pip or setuptools at this time).  Then:

```sh
% python setup.py install
````

Running
-------
Basic usage is as follows:

```sh
% quiver aligned_reads{.cmp.h5, .bam, .fofn, or .xml}    \
>     -r reference{.fasta or .xml} -o variants.gff       \
>     -o consensus.fasta -o consensus.fastq
```

in this example we perform haploid consensus and variant calling on
the mapped reads in the ``aligned_reads.bam`` which was aligned to
``reference.fasta``.  The ``reference.fasta`` is only used for
designating variant calls, not for computing the consensus.  The
consensus quality score for every position can be found in the output
FASTQ file.


Documentation
-------------

- [More detailed installation instructions](https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/HowToQuiver.rst)
- [Quiver FAQ](https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/QuiverFAQ.rst)
- [variants.gff spec](https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/VariantsGffSpecification.rst)
- [CHANGELOG](https://github.com/PacificBiosciences/GenomicConsensus/blob/master/CHANGELOG)
