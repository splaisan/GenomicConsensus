GenomicConsensus (quiver, arrow) [![Circle CI](https://circleci.com/gh/PacificBiosciences/GenomicConsensus.svg?style=svg)](https://circleci.com/gh/PacificBiosciences/GenomicConsensus)
-------------------------

The ``GenomicConsensus`` package provides the ``variantCaller`` tool,
which allows you to apply the Quiver or Arrow algorithm to the mapped
PacBio reads to get consensus and variant calls.

Background on Quiver and Arrow
------------------------------

*Quiver* is the legacy consensus model based on a conditional random field approach.  Over the
years it has proven difficult to train and develop, so we are phasing it out in favor of the
new model, Arrow.

*Arrow* is an improved consensus model based on a more straightforward hidden Markov model
approach.

Quiver is supported for PacBio RS data.  Arrow is supported for PacBio Sequel data
and RS data with P6-C4 chemistry.


Getting ``GenomicConsensus``
----------------------------
Casual users should get ``GenomicConsensus`` from the [SMRTanalysis software bundle](http://www.pacb.com/support/software-downloads/).

Instructions for those who want to build the software themselves are [here](./doc/Installation.rst).

Running
-------
Basic usage is as follows:

```sh
% quiver aligned_reads{.cmp.h5, .bam, .fofn, or .xml}    \
>     -r reference{.fasta or .xml} -o variants.gff       \
>     -o consensus.fasta -o consensus.fastq
```

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


More documentation
------------------

- [Detailed installation instructions](./doc/Installation.rst)
- [FAQ](./doc/FAQ.rst)
- [variants.gff spec](./doc/VariantsGffSpecification.rst)
- [CHANGELOG](./CHANGELOG)
