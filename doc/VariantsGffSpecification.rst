
``variants.gff`` File Format (Version 1.4)
============================================

As of this version, ``variants.gff`` is our primary variant call file
format.  The ``variants.gff`` file is based on the `GFFv3 standard`_.
The GFFv3 standard describes a tab-delimited plain-text file
meta-format for describing genomic "features."  Each gff file consists
of some initial "header" lines supplying metadata, and then a number
of "feature" lines providing information about each identified
variant.

The GFF Coordinate System
-------------------------

All coordinates in GFF files are 1-based, and all intervals ``start,
end`` are understood as including both endpoints.

Headers
-------

The ``variants.gff`` file begins with a block of metadata headers,
which looks like the following:

::

    ##gff-version 3
    ##pacbio-variant-version 1.4
    ##date Tue Feb 28 17:44:18 2012
    ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
    ##source GenomicConsensus v0.1.0
    ##source-commandline callVariants.py --algorithm=plurality aligned_reads.cmp.h5 -o variants.gff
    ##sequence-header ref000022 EGFR_Exon_23
    ##sequence-region ref000022 1 189
    ##sequence-region ref000022 200 235
    ##sequence-header ref000023 EGFR_Exon_24
    ##sequence-region ref000023 1 200

The ``source`` and ``source-commandline`` describe the name and
version of the software generating the file.  ``sequence-header``
describes a mapping between the full name of a reference group (i.e. a
reference contig), and the local identifier that will be used to
denote that reference group in the file.  Here, the full name is
``EGFR_Exon_23``, while the local identifier is ``ref00022``.  The
local identifier will always be purely alphanumeric.

There can be many reference groups and regions present in the file.
They will always represented as the last portion of the header block,
in the order *sequence-header*, corresponding *sequence-regions*.

``pacbio-variant-version`` reflects the specification version that the
file contents should adhere to.



Feature lines
-------------

After the headers, each line in the file describes a genomic
*feature*; in this file, all the features are potential variants
flagged by the variant caller.  The general format of a variant line
is a 9-column (tab-delimited) record, where the first 8 columns
correspond to fixed, predefined entities in the GFF standard, while
the 9th column is a flexible semicolon-delimited list of mappings
``key=value``.

The 8 predefined columns are as follows:

+------+-------+--------------------------------+------------------+
|Column|Name   |Description                     |Example           |
|Number|       |                                |                  |
+------+-------+--------------------------------+------------------+
|1     |seqId  |The full FASTA header for the   |``lambda_NEB3011``|
|      |       |reference contig.               |                  |
|      |       |                                |                  |
+------+-------+--------------------------------+------------------+
|2     |source |(unused; always populated with  |``.``             |
|      |       |``.``)                          |                  |
+------+-------+--------------------------------+------------------+
|3     |type   |the type of variant.  One of    |``substitution``  |
|      |       |``insertion``, ``deletion``, or |                  |
|      |       |``substitution``.               |                  |
|      |       |                                |                  |
+------+-------+--------------------------------+------------------+
|4     |start  |1-based start coordinate for the|200               |
|      |       |variant.                        |                  |
+------+-------+--------------------------------+------------------+
|5     |end    |1-based end coordinate for the  |215               |
|      |       |variant.  start<=end always     |                  |
|      |       |obtains, regardless of strand.  |                  |
+------+-------+--------------------------------+------------------+
|6     |score  |unused; populated with ``.``    |``.``             |
+------+-------+--------------------------------+------------------+
|7     |strand |unused; populated with ``.``    |``.``             |
|      |       |                                |                  |
+------+-------+--------------------------------+------------------+
|8     |phase  |unused; populated with ``.``    |``.``             |
+------+-------+--------------------------------+------------------+


The attributes in the 9th (final) column are as follows:

+--------------+----------------------------+-----------------+
|Key           |Description                 |Example          |
|              |                            |value            |
+--------------+----------------------------+-----------------+
|``length``    |the length of the variant   |``1``            |
|              |                            |                 |
|              |                            |                 |
|              |                            |                 |
+--------------+----------------------------+-----------------+
|``coverage``  |the read coverage of the    |``42``           |
|              |variant site (not the       |                 |
|              |variant itself)             |                 |
+--------------+----------------------------+-----------------+
|``confidence``|the phred-scaled probability|``37``           |
|              |that the variant is real,   |                 |
|              |rounded to the nearest      |                 |
|              |integer and truncated at 93 |                 |
+--------------+----------------------------+-----------------+
|``reference`` |the reference base or bases |``T``            |
|              |for the variant site        |                 |
+--------------+----------------------------+-----------------+
|``variantSeq``|the read base or bases      |``T``            |
|              |corresponding to the variant|(haploid);       |
|              |                            |``T,C``          |
|              |                            |(diploid)        |
+--------------+----------------------------+-----------------+
|``zygosity``  |One of ``heterozygous``, or | ``homozygous``  |
|              |``homozygous``, to denote   |                 |
|              |the zygosity of the variant.|                 |
|              | This attribute is not used |                 |
|              |for haploid analyses.       |                 |
|              |                            |                 |
+--------------+----------------------------+-----------------+
|``frequency`` |the read coverage of the    |``13``           |
|              |variant itself              |                 |
+--------------+----------------------------+-----------------+


The attributes may be present in any order.

The four types of variant we support are as follows. *(Recall that the
field separator is a tab, not a space.)*

1. Insertion.  Example::

    ref00001 . insertion 8 8 . . . length=1;variantSeq=G;confidence=22;coverage=18

  For insertions, start==end, and the insertion event is understood as
  taking place *following* the reference position `start`.  `length`
  corresponds to the length of the insertion.  The `reference` attribute
  is omitted for insertions.

2. Deletion.  Example::

    ref00001 . deletion 348 349 . . . length=1;reference=G;confidence=39;coverage=25

  For deletions, the length is the length of reference sequence
  deleted.  The `variantSeq` attribute is omitted for deletions.


3. Substitution.  Example::

    ref000001 . substitution 100 102 . . . length=3;reference=GGG;variantSeq=CCC;confidence=50;coverage=20


Compression
-----------

The gff metaformat is verbose, so for practical purposes we will gzip
encode ``variants.gff`` files as ``variants.gff.gz``.  Consumers of
the variant file should be able to read it in either form.


Other file formats
------------------

The VCF and BED standards describe variant-call specific file formats.
We can currently translate `variants.gff` files to these formats, but
they are not the primary output of the variant callers.


.. _GFFv3 standard: http://www.sequenceontology.org/gff3.shtml
