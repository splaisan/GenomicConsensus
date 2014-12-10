Variant Caller Validation Specification
=======================================

Created: 2012/02/20
Author: jdrake

Synopsis
--------
There are several algorithms implemented for detecting SNPs in the data using alignments
against a known reference.  These include Evicons (PacBio), GATK (Broad) and, most recently,
GenomicConsensus (PacBio).  The first was built back in the days of 70% accurate reads
and is quickly becoming deprecated, though currently used only as our haploid caller (though it
does have diploid calling functionality).  The second is part of a comprehensive tool kit
that provides better diploid calling and is currently used as such in the secondary pipeline.
The third is the most recent incarnation and the heir apparent going forward.

There are no metrics built to measure, for example, the sensitivy and specificity of these
algorithms, and thus are difficult to evaluate against eachother.  Ostensibly, since we're
closer to the data, we should be able to better tune an algorithm to maximize true +/-
variant calls.  This exercise will create datasets to generate ROC curves and potentially
other user metrics to properly evaluate algorithms.

Workflow
--------
A dataset will be generated using a set of mutated reference genomes to align real reads
against.  Each mutated reference will have a list of 'ground-truth' mutations associated
with them.  The alignments will then be processed by each of the candidate variant caller
algorithms and their results evaluated against the ground truth for true +/-.

1. Generate mutated reference(s)
2. Align reads to each mutated reference
3. Run variant callers using alignments
4. Evaluate calls vs ground truth
5. Generate metrics
6. Repeat steps 3 - 5 as necessary

*NOTE*: The mutated references could be generated on the fly using the mutation information,
thus, obviating the need to save all the mutated references.

The automation of the workflow should eventually be packaged and deployed into the Millhouse
framework. Initially, it can configured to run from the command line on the secondary cluster
and the smrtpipe infrastructure.

Mutation Sets
-------------
Starting with lambda (well represented amongst currently available runs), generate a set, M, of
mutated lambda references with n randomly generated point mutations within each m mutated genome.
Each point mutation p will be one of P = {Substitution(s), Insertion(i), Deletion(d)}. Locations will
be associated with each mutation as a 1-based offset using the wild-type genome (w) coordinates.

Offsets are stored in 0-based coordinates, but displayed and manipulated using a 1-based coordinate
system because that's what GFF, SAM and Tablet uses.

Mutation sets will be stored in a file per reference mutated.  The file will be used as input to
mutate genomes on the fly just prior to alignment as well as input to the validation procedure.
GFF files could be generated from them fairly easily if, though probably not, necessary.  Multiple
versions of this file could co-exist for the same reference.

The format of the mutations file is extremely simple and compact making it suitable for source
control.  For simplicity, we'll use the python pickle protocol vers 2.

mutation = {
    offs, # offset in the wild-type genome (1-based)
    typ,  # type of mutation
    wild, # base(s), wild-type strain
    mut   # base(s), the mutated strain
}

Comparison
~~~~~~~~~~

Alignments play a big role in this.  Homopolymer regions are treated slightly differently when
QV's are involved.  Without them, affine gaps are used which push gaps to the left, e.g.,

GAATGAAGCCG
GAATG-AGCCG

These degenerate cases will be handled by collapsing the homopolymer aligment gaps to the left before
comparing.  See BLASR subsection for more detail.

Some alignments against a substitution generate this type of call (some field removed for brevity):
deletion 10201 length=1;confidence=22;coverage=16;reference=G
insertion 10201 length=1;variantSeq=T;confidence=23;coverage=16

This happens when the substitution forms a homopolymer.  These will be collapsed and labelled substitutions
during comparisons.

Does it use CCS reads by default, QV values?  The production version does not use CCS reads but does use
QV values.


Data Sets
---------
Key things to pay attention to in datasets:
- coverage level
- quality
- location of mutation (e.g., homopolymer)
- nature of mutation (e.g., insertion followed by deletion)

Start with a positive control using an unmutated reference.  Zero mutations should be found.

Metrics
-------
Confusion matrix
ROC Curves (using quality scores)
QQ Plot

Notes
-----
http://www.sequenceontology.org/gff3.shtml
http://smrtwiki/wiki/SMRTpipe
http://web/~mhsieh/dag/index.html

Evicons
~~~~~~~
Top level module src: //depot/software/assembly/pbpy/pbpy/smrtpipe/modules
Evicons smrtpipe modules: P_Consensus, P_ConsensusAlgorithm
    wraps runChunkConnectPost.py (same dir) which ...
    wraps eviConsWrapper.py (../../../bin/) which ...
    wraps jConnectPost[.sh] (../../../../seymour/dist2/analysis/bin/) which ...
    wraps a call an evicons jar file
*Un*-wrapping this may be more cumbersome than generating the appropriate inputs to the module.

GATK
~~~~
P_GATKVC

Uses the UnifiedGenotyper, TableRecalibration, CountCovariates components

Uses BAM inputs, generated after alignment (blasr)
//depot/software/bioinformatics/tools/pbsamtools

BLASR
~~~~~
Running blasr to get a cmp.h5 file (super basic, with crappy alignments)::

> compareSequences.py --algorithm=blasr --h5fn=aligned.cmp.h5 input.fofn refdir

More productiony way::

> compareSequences.py --info --useGuidedAlign --algorithm=blasr --nproc=6  --noXML --h5mode=w \
    --h5fn=control.cmp.h5 --minAccuracy=0.75 --minLength=50  -x -minMatch 12 -x -bestn 1 -x -minPctIdentity 70.0 \
    --regionTable=/mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/jobs/037/037285/data/filtered_regions.chunk001of002.fofn \
    input.fofn /mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/references/lambda


`refdir` is a directory containing a set of information related to a reference sequence.
The key files appear to be <reference>.fa and reference.info.xml.  It can work with just
a fasta file, but will produce a cmp.h5 that breaks evicons (reference length is 0). There
is a utility to generate these ref dirs:

/mnt/secondary/Smrtpipe/builds/Assembly_Mainline_Nightly_LastSuccessfulBuild/analysis/bin/referenceUploader


Validation tests could be source controlled under the siv tree, given they're likely to
transition into that group eventually (//depot/software/assembly/siv-test/...)

Using what we've already got:
//depot/software/assembly/siv-test/module-test/bin/
- mutateRef.py (?)
- evalVariantCalls.py (?)

