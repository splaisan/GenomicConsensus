1.3.3 Enhancements
==================

Bug 20100
---------
**Genomic Consensus to support rare variant calling**

Adds the ability to to call rare variants.  Rare variants are defined here as
mutations detected at a 1% < frequency < 50%.  There is an initial minimum
coverage requirement set at 500x.  The information provided by this feature
will be limited to deviations from the reference.

This will limited to SNPs only.  Indels will be ignored.

Codon-aware filtering could easily be applied to the output as a post-processing
step.  It could also be used *in situ* as an additional filtering mechanism to
reduce potentially noisy output.  We can start with the post-processing option
(easy) then evolve towards being codon aware as necessary.

This functionality will be optional.

*Inputs*:

    A column-oriented set of short (~1 - 3 bases) sequence calls and their
    corresponding frequencies.

*Outputs*:

    A GFF-style record per rare variant call. See VariantsGffSpecification for
    standard format. This feature will augment the standard record with a new
    key/value pair indicating frequency. Example: freq=0.10.  There may be more
    than one variant per reference position.

    Please note that no consensus file(s) will be generated for rare variants,
    though enough information is provided to build one in a separate
    tools/module.

Bug 20628
---------
**Add support for detecting and reporting correlated mutations**

Provides support for determing whether or not sets of mutations are co-located.
This only includes SNPs, not indels.  Correlations may only be established
using overrlapping reads, i.e., gaps not supportable. Correlations will have
some confidence metric associated with them (TBD).  This functionality may also
be combined with rare variant calling output.

The guiding use-case for this feature is the BCR-ABL project, which targeted an
863 bp region of the human genome (11,000x coverage). `BCR-ABL Project Details`_.

This functionality will be optional.

*Inputs*:

    CCS-based variant calls at each position including read information: ID,
    start, stop.  'start' and 'stop' are in (+) strand genomic coordinates.

*Outputs*:

    A table (possibly) of correlated mutations that could look like:

    =====  =======  =====  ===================
    Freq   # Reads  Conf   Mutations
    =====  =======  =====  ===================
    40.4%  4,321    40     123a, 140c
    30.3%  3,210    30     50t, 350a
    20.2%  2,500    20     1400g, 1500c, 1550t
    =====  =======  =====  ===================

    We may also choose to include an output of read IDs associated with reported
    sets of co-located mutations.  Formats TBD.

.. _BCR-ABL Project Details: http://web/~mbrown/workspace2011Q4/bcrAblASHRuns/
