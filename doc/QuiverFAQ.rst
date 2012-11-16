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
For example, consider this example::

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
additional /QV/ covariates that the basecaller provides.  Using these
covariates provides additional information tat

Quiver does not use the alignment provided by the mapper (BLASR,
typically), except for determining how to group reads together at the
gross level.  It implicitly performs its own realignment, so it is
highly sensitive to all variant types, including indels---for example
it resolves the example above with ease.

The name **Quiver** reflects a consensus-calling algorithm that is
`QV-aware`.

How do I run Quiver?
--------------------
For general instructions on installing and running, see the
HowToQuiver_ document.



What does Quiver put in its output files?
-----------------------------------------
There are three output files that you can get from Quiver:

    - A consensus *FASTA* file containing the consensus sequence
    - A consensus *FASTQ* file containing the consensus sequence with quality annotations
    - A variants *GFF* file containing a filtered, annotated list of variants identified

It is important to note that the variants included in the output
variants GFF file are /filtered/ by coverage and quality, so not all
variants that are apparent in comparing the reference to the consensus
FASTA output will correspond to variants in the output variants GFF
file.

To enable all output files, you can run, for example::

     % quiver -j16 aligned_reads.cmp.h5 -r ref.fa \
         -o consensus.fa                          \
         -o consensus.fq                          \
         -o variants.gff


The extension is used to determine the output file format.


What does it mean that Quiver's consensus is /de novo/?
-------------------------------------------------------
It is /de novo/ in the sense that, in its default configuration, the
reference and the reference alignment are not used to inform the
consensus output.  Only the reads factor into the determination of the
consensus.


What is Quiver's accuracy?
--------------------------
Quiver's expected accuracy is a function of coverage.  Using the C2
chemistry, the accuracy we typically achieve is >Q50@60x; >Q40@40x;
>Q30@20x; >Q20@5x.  The "Q" values we refer to are Phred-scaled
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



Why is Quiver making errors in some region?
-------------------------------------------
The most likely cause for /true/ errors made by Quiver is that the
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
ensure that the region was called at high confidence; then see



What is `MapQV` and why is it important?
----------------------------------------
`MapQV` is a single scalar Phred-scaled QV per aligned read, that
reflects the mapper's degree of certainty that the read aligned to
/this/ part of the reference and not some other.  Unambigously mapped
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
these regions, Quiver will, by default, opt to "no-call" these
regions---putting "N" in consensus output there.

You have two options when confronted with this problem.

First, you can turn off the `MapQV` filter entirely.  In this case,
the consensus for each instance of a genomic repeat will be calculated
using reads that may actually be from other instances of the repeat,
so the exact trustworthiness of the consensus in that region may be
suspect.

Second, if you believe that your original reference is relatively
accurate, you can have Quiver insert the reference bases into the
consensus output in these regions of low effective coverage.


How can I turn off the `MapQV` filter and why would I want to?
--------------------------------------------------------------
You can disable the `MapQV` filter using the flag `--mapQvThreshold=0`
(shorthand: `-m 0`).  You might want to do this in de novo assembly
projects, but it is not recommended for variant calling applications.


How can I make Quiver output the original reference in areas of low coverage?
-----------------------------------------------------------------------------
When Quiver is confronted with a region where effective coverage is so
low that high-quality consensus cannot be produced, it  has two options:

  - "no-call" the region, filling the consensus with "N" bases and not
    making any variant calls in the region.  This can be enabled by
    `--noEvidenceConsensusCall=nocall`---but it is the default.  This
    is suitable when you believe the original reference may be
    inaccurate.

  - transfer the reference sequence in the region into the consensus
    output, and call no variants in the region.  This can be enabled
    by `--noEvidenceConsensusCall=reference`.  This would be a good
    option if you bleieve your reference to be relatively accurate.


How do I inspect or validate the variant calls made by Quiver?
--------------------------------------------------------------
When in doubt, it is easiest to inspect the region in a tool like
SMRTView®, which enables you to view the reads aligned to the region.
Deletions and substitutions should be fairly easy to spot; to view
insertions, right-click on the reference base and select "View
Insertions Before...".

Another approach is to use the `--dumpEvidence` flag, which will
output a directory for each window surrounding a called variant,
containing all the reads clipped to the window.  You can use an
independent consensus calling approach or build and view a multiple
alignment from these reads.


What are the filtering parameters that Quiver uses?
---------------------------------------------------

Quiver limits read coverage, filters reads by `MapQV`, and filters
variants by quality and coverage.

 - The overall read coverage used to call consensus in every window is
   100x by default, but can be changed using ``-Xvalue``.
 - The `MapQV` filter, by default, removes reads with MapQV < 20.  This
   is configured using ``--mapQvThreshold=value`` / ``-m value``
 - Variants are only called if the read coverage of the site exceeds
   11x, by default---this is configurable using ``-x value``.
   Further, they will not be called if the confidence (Phred-scaled)
   does not exceed 20---configurable using ``-q value``.



What is the best way to call consensus on an amplicon dataset?
--------------------------------------------------------------
In an amplicon dataset, focused regions of a genome have been
amplified, ideally with minimal off-target amplification.  If you
provide Quiver a reference that is the full genome, not just the
amplified regions, it will get tripped up by the large regions

To avoid this problem, it is best to split out each amplicon region of
the reference into its own reference contig.



What happens when my sample is a mixture, or diploid?
-----------------------------------------------------
At present, Quiver assumes a haploid sample, and the behavior of *Quiver* on sample mixtures or
diploid/polyploid samples is /undefined/.  The program will not crash,
but the output results are not guaranteed to accord with any one of
the haplotypes in the sample, as opposed to a potential patchwork.  We
are working on improvements for the 2.0 release.


.. _HowToQuiver: https://github.com/PacificBiosciences/GenomicConsensus/blob/master/doc/HowToQuiver.rst
