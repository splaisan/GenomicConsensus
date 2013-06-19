Quiver FAQ
==========

EviCons? GenomicConsensus? Quiver? Plurality?  I'm confused!
------------------------------------------------------------
**GenomicConsensus** is the current PacBio consensus and variant
calling suite.  It contains a main program, ``variantCaller.py``,
which provides two consensus / variant calling algorithms: **Plurality**
and **Quiver**.  You can run these algorithms by calling
``variantCaller.py --algorithm=[quiver|plurality]`` or through the
convenience wrapper scripes ``quiver`` and ``plurality``.

**EviCons** was the previous generation PacBio variant caller and was
deprecated in release 1.3.1.

Not to confuse you, but there is also a separate package called
**ConsensusCore**.  This is a package containing a C++ library where all
the computation behind Quiver is done.  You shouldn't have to worry
about it once it's installed.


What is Plurality?
------------------
**Plurality** is a very simple variant calling algorithm: stack up the
aligned reads (alignment as produced by BLASR, or alternate mapping
tool), and for each column under a reference base, call the most
abundant (i.e., the plurality) read base (or bases, or deletion) as
the consensus at that reference position.


Why is Plurality a weak algorithm?
----------------------------------
Plurality does not perform any local realignment.  This means it is
heavily biased by the alignment produced by the mapper (BLASR,
typically).  It also means that it is insensitive at detecting indels.
Consider this example::

    Reference    AAAA
                 ----
      Aligned    A-AA
        reads    AA-A
                 -AAA
                 ----
    Plurality    AAAA
    consensus

Note here that every read has a deletion and the correct consensus
call would be "AAA", but due to the mapper's freedom in gap-placement
at the single-read level, the plurality sequence is "AAAA"---the
deletion is missed.  Local realignment, which plurality does not do,
but which could be considered as implicit in the Quiver algorithm,
essentially pushes the gaps here to the same column, thus identifying
the deletion.

What is Quiver?
---------------
**Quiver** is a more sophisticated algorithm that finds the maximum
likelihood template sequence given PacBio reads of the template.  We
model PacBio reads using a conditional random field approach that
prescribes a probability to a read given a template sequence.  In
addition to the base sequence of each read, Quiver uses several
additional *QV* covariates that the basecaller provides.  Using these
covariates provides additional information about each read, allowing
more accurate consensus calls.

Quiver does not use the alignment provided by the mapper (BLASR,
typically), except for determining how to group reads together at the
gross level.  It implicitly performs its own realignment, so it is
highly sensitive to all variant types, including indels---for example
it resolves the example above with ease.

The name **Quiver** reflects a consensus-calling algorithm that is
`QV-aware`.

How do I run Quiver?
--------------------

If you have
For general instructions on installing and running, see the
HowToQuiver_ document.



What does Quiver put in its output files?
-----------------------------------------
There are three output files that you can get from Quiver:

- A consensus *FASTA* file containing the consensus sequence
- A consensus *FASTQ* file containing the consensus sequence with quality annotations
- A variants *GFF* file containing a filtered, annotated list of variants identified

It is important to note that the variants included in the output
variants GFF file are *filtered* by coverage and quality, so not all
variants that are apparent in comparing the reference to the consensus
FASTA output will correspond to variants in the output variants GFF
file.

To enable all output files, you can run, for example::

     % quiver -j16 aligned_reads.cmp.h5 -r ref.fa \
         -o consensus.fa                          \
         -o consensus.fq                          \
         -o variants.gff


The extension is used to determine the output file format.


What does it mean that Quiver's consensus is *de novo*?
-------------------------------------------------------
It is *de novo* in the sense that the reference and the reference
alignment are not used to inform the consensus output.  Only the reads
factor into the determination of the consensus.

The only time the reference sequence is used to make consensus calls,
when the ``--noEvidenceConsensusCall`` flag is set to ``reference`` or
``lowercasereference`` (the default), is when there is no effective
coverage in a genomic window, so Quiver has no evidence for computing
consensus.  The purist can set ``--noEvidenceConsensusCall=nocall`` to
avoid using the reference even in zero coverage regions.


What is Quiver's accuracy?
--------------------------
Quiver's expected accuracy is a function of coverage.  Using the C2
and P4-C2 chemistries, our nominal expected accuracy levels are as
follows:

+--------+---------+
|Coverage|Expected |
|        |consensus|
|        |accuracy |
+========+=========+
|10x     | > Q30   |
+--------+---------+
|20x     | > Q40   |
+--------+---------+
|40x     | > Q50   |
+--------+---------+
|60-80x  | ~ Q60   |
+--------+---------+

The "Q" values we refer to are Phred-scaled
quality values::

   q = -10 log_10 p_error

so for instance Q50 corresponds to a p_error of 0.00001---an accuracy
of 99.999%.

We are working to lower the coverage requirements needed to achieve
high accuracy.


What are these QVs that Quiver uses?
------------------------------------
Quiver uses additional QV tracks provided by the basecaller.  I like
to think of these QVs as little breadcrumbs that are left behind by
the basecaller to help identify positions where it was likely that
errors of given type occurred.  Formally, the QVs for a given read are
vectors of the same length as the number of bases called; the QVs we
use are as follows:

  - DeletionQV
  - InsertionQV
  - MergeQV
  - SubstitutionQV
  - DeletionTag

To find out if your cmp.h5 file is loaded with these QV tracks, run the command
::

    % h5ls -rv aligned_reads.cmp.h5

and look for the QV track names in the output.  If your cmp.h5 file is
lacking some of these tracks, Quiver will still run, though it will
issue a warning that its performance will be suboptimal.


Why is Quiver making errors in some region?
-------------------------------------------
The most likely cause for *true* errors made by Quiver is that the
coverage in the region was low.  If you only have 5x coverage over a
1000-base region, you would expect 10 errors in that region.

It is important to understand that the effective coverage available to
Quiver is not the full coverage apparent in plots---Quiver and
Plurality both filter out ambiguously mapped reads by default.  The
remaining coverage after filtering is called the /effective coverage/.
See the next section for discussion of `MapQV`.

If you have verified that there is high effective coverage in region
in question, it is highly possible---given the high accuracy Quiver
can achieve---that the apparent errors you are observing actually
reflect true sequence variants.  Inspect the FASTQ output file to
ensure that the region was called at high confidence; if an erroneous
sequence variant is being called at high confidence, please report a
bug to us.


What does Quiver do for genomic regions with no effective coverage?
-------------------------------------------------------------------
For regions with no effective coverage, no variants are outputted, and
the FASTQ confidence is 0.

The output in the FASTA and FASTQ consensus sequence tracks is
dependent on the setting of the ``--noEvidenceConsensusCall`` flag.
Assuming the reference in the window is "ACGT", the options are:

+---------------------------------------------+---------+
|``--noEvidenceConsensusCall=...``            |Consensus|
|                                             |output   |
+=============================================+=========+
|``nocall`` (default in 1.4)                  |NNNN     |
+---------------------------------------------+---------+
|``reference``                                |ACGT     |
+---------------------------------------------+---------+
|``lowercasereference`` (new post 1.4, and the|         |
|default)                                     |acgt     |
+---------------------------------------------+---------+




What is `MapQV` and why is it important?
----------------------------------------
`MapQV` is a single scalar Phred-scaled QV per aligned read, that
reflects the mapper's degree of certainty that the read aligned to
*this* part of the reference and not some other.  Unambigously mapped
reads will have a high `MapQV` (typically 255), while a read that was
equally likely to have come from two parts of the reference would have
a `MapQV` of 3.

`MapQV` is pretty important when you want highly accurate variant
calls.  Quiver and Plurality both filter out aligned reads with a
MapQV below 20 (by default), so as not to call a variant using data of
uncertain genomic origin.

This can cause problems when you are using Quiver to get a consensus
sequence.  If your genome contains long (relative to your library
insert size) highly-similar repeats, the effective coverage (after
`MapQV` filtering) may be reduced in the repeat regions---we term
these `MapQV` dropouts.  If the coverage is sufficiently reduced in
these regions, Quiver will not call consensus in these regions---see
`What does Quiver do for genomic regions with no effective coverage?`_.

If you want to use ambiguously mapped reads in computing a consensus
for a denovo assembly, you can turn off the `MapQV` filter entirely.
In this case, the consensus for each instance of a genomic repeat will
be calculated using reads that may actually be from other instances of
the repeat, so the exact trustworthiness of the consensus in that
region may be suspect.  The next section describes how to disable the
`MapQV` filter.


How can I turn off the `MapQV` filter and why would I want to?
--------------------------------------------------------------
You can disable the `MapQV` filter using the flag
``--mapQvThreshold=0`` (shorthand: ``-m=0``).  If you are running your
Quiver job via SMRTportal, this can be done by unchecking the "Use
only unambiguously mapped reads" option. You might want to do this in
de novo assembly projects, but it is not recommended for variant
calling applications.


How do I inspect or validate the variant calls made by Quiver?
--------------------------------------------------------------
When in doubt, it is easiest to inspect the region in a tool like
SMRTView®, which enables you to view the reads aligned to the region.
Deletions and substitutions should be fairly easy to spot; to view
insertions, right-click on the reference base and select "View
Insertions Before...".


What are the filtering parameters that Quiver uses?
---------------------------------------------------

Quiver limits read coverage, filters reads by `MapQV`, and filters
variants by quality and coverage.

- The overall read coverage used to call consensus in every window is
  100x by default, but can be changed using ``-X=value``.
- The `MapQV` filter, by default, removes reads with MapQV < 20.  This
  is configured using ``--mapQvThreshold=value`` / ``-m=value``
- Variants are only called if the read coverage of the site exceeds
  5x, by default---this is configurable using ``-x=value``.
  Further, they will not be called if the confidence (Phred-scaled)
  does not exceed 40---configurable using ``-q=value``.


What happens when my sample is a mixture, or diploid?
-----------------------------------------------------
At present, Quiver assumes a haploid sample, and the behavior of
*Quiver* on sample mixtures or diploid/polyploid samples is
*undefined*.  The program will not crash, but the output results are
not guaranteed to accord with any one of the haplotypes in the sample,
as opposed to a potential patchwork.  We are working on improvements
for the 2.0 release.


.. _HowToQuiver: https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/HowToQuiver.rst
