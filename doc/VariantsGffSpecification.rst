
``variants.gff`` File Format (Version 2.1)
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
    ##pacbio-variant-version 2.1
    ##date Tue Feb 28 17:44:18 2012
    ##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12
    ##source GenomicConsensus v0.1.0
    ##source-commandline callVariants.py --algorithm=plurality aligned_reads.cmp.h5 -r spinach.fasta -o variants.gff
    ##source-alignment-file /home/popeye/data/aligned_reads.cmp.h5
    ##source-reference-file /home/popeye/data/spinach.fasta
    ##sequence-region EGFR_Exon_23 1 189
    ##sequence-header EGFR_Exon_24 1 200

The ``source`` and ``source-commandline`` describe the name and
version of the software generating the file.
``pacbio-variant-version`` reflects the specification version that the
file contents should adhere to.

  The ``sequence-region`` headers describe the names and extents of
the reference groups (i.e. reference contigs) that will be refered to
in the file.  The names are the same as the full FASTA header.

``source-alignment-file`` and ``source-reference-file`` record
absolute paths to the primary input files.


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
|``coverage``  |the read coverage of the    |``42``           |
|              |variant site (not the       |                 |
|              |variant itself)             |                 |
+--------------+----------------------------+-----------------+
|``confidence``|the phred-scaled probability|``37``           |
|              |that the variant is real,   |                 |
|              |rounded to the nearest      |                 |
|              |integer and truncated at 93 |                 |
+--------------+----------------------------+-----------------+
|``reference`` |the reference base or bases |``T``, ``.``     |
|              |for the variant site.  May  |                 |
|              |be ``.`` to represent a     |                 |
|              |zero-length substring (for  |                 |
|              |insertion events)           |                 |
+--------------+----------------------------+-----------------+
|``variantSeq``|the read base or bases      |``T``            |
|              |corresponding to the        | (haploid);      |
|              |variant. ``.`` encodes a    |``T/C``, ``T/.`` |
|              |zer-length string, as for a | (heterozygous)  |
|              |deletion.                   |                 |
+--------------+----------------------------+-----------------+
|``frequency`` |the read coverage of the    |``13``           |
|              |variant itself; for         | (haploid)       |
|              |heterozygous variants, the  |                 |
|              |frequency of both observed  |``15/12``        |
|              |alleles.  This is an        | (heterozygous)  |
|              |optional field.             |                 |
+--------------+----------------------------+-----------------+


The attributes may be present in any order.

The three types of variant we support are as follows. *(Recall that the
field separator is a tab, not a space.)*

1. Insertion.  Examples::

    ref00001 . insertion 8 8 . . . reference=.;variantSeq=G;confidence=22;coverage=18;frequency=10
    ref00001 . insertion 19 19 . . . reference=.;variantSeq=G/.;confidence=22;coverage=18;frequency=7/5

  For insertions, start==end, and the insertion event is understood as
  taking place *following* the reference position `start`.

2. Deletion.  Examples::

    ref00001 . deletion 348 349 . . . reference=G;variantSeq=.;confidence=39;coverage=25;frequency=20
    ref00001 . deletion 441 443 . . . reference=GG;variantSeq=GG/.;confidence=39;coverage=25;frequency=8/8

3. Substitution.  Examples::

    ref000001 . substitution 100 102 . . . reference=GGG;variantSeq=CCC;confidence=50;coverage=20;frequency=16
    ref000001 . substitution 200 201 . . . reference=G;variantSeq=G/C;confidence=50;coverage=20;frequency=10/6



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
