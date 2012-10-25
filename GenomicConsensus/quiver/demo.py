## Let's call consensus using Quiver.

# To run this script, you will need:
#
#    - python2.7
#    - a virtualenv containing current (mainline) versions of
#      * GenomicConsensus
#      * ConsensusCore
#      * pbcore
#
# If you just want to see what the output looks like, see
# [here](example_output.txt).


# Load some bits we'll need.
import random, numpy as np, sys
from collections import Counter
from pbcore.io import CmpH5Reader, FastaReader

import ConsensusCore as cc
from GenomicConsensus.quiver.utils import *
from GenomicConsensus.quiver.model import *


# Set up some parameters.  Don't worry about them too much for now.
# We only use 11x POA coverage as the POA performance seems to degrade
# beyond that.
PARAMS          = AllQVsModel.trainedParams1().qvModelParams
POA_COVERAGE    = 11
QUIVER_COVERAGE = 30

# Some lisp functions we want
fst = lambda t: t[0]
snd = lambda t: t[1]

# Quiver requires QVs that aren't propagated to the cmp.h5 files by
# default in C2.  You can get them by running the Resequencing_QVs
# workflow in smrtportal.  Here are a couple canned examples.
if len(sys.argv) > 1 and sys.argv[1] == "ecoli":
    refName   = "E. coli"
    refFasta  = "/mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/references/ecoli_mutated/sequence/ecoli_mutated.fasta"
    inputFile = "/mnt/secondary/Smrtanalysis/userdata/jobs/044/044601/data/aligned_reads.cmp.h5"
else:
    refName   = "Lambda"
    refFasta  = "/mnt/secondary/Smrtanalysis/opt/smrtanalysis/common/references/lambda/sequence/lambda.fasta"
    inputFile = "/mnt/secondary/Smrtanalysis/userdata/jobs/038/038537/data/aligned_reads.cmp.h5"

c = CmpH5Reader(inputFile)
refString = next(iter(FastaReader(refFasta))).sequence

# This is the heart of the matter.  For demonstration purposes, we're
# processing the genome in a single loop in a single thread
# (/process).
#
# What we do here is: loop over 1000-base windows of the reference,
# taking all reads that span each reference window and clipping them
# to the bounds.  Note that for this demo we are not using any reads
# that do not completely span a window.
table = Counter()

for start in xrange(0, len(refString), 1000):
    stop = min(start + 1000, len(refString))
    win = (start, stop)
    print "%s [%d, %d)" % (refName, start, stop)

    reads = c.readsInRange(1, *win)
    spanningReads = [a for a in reads
                     if (a.referenceStart <= start and stop <= a.referenceEnd)]
    clippedReads = [a.clippedTo(*win) for a in spanningReads]

    # Bail out if there are no spanning reads
    if clippedReads == []:
        print "Inadequate spanning coverage"
        table += Counter({"Skipped" : 1000})
        continue

    # We calculate a quick consensus estimate using the
    # Partial Order Aligner (POA), which typically gets us to > 99%
    # accuracy.
    strs = [a.read(orientation="genomic", aligned=False) for a in clippedReads]
    poa = cc.PoaConsensus.FindConsensus(strs[:POA_COVERAGE])
    ga = cc.Align(refString[start:stop], poa.Sequence())
    print "  POA:      %0.3f  (%d errors)" % (ga.Accuracy(), ga.Errors())

    roughConsensus = poa.Sequence()

    # Load the reads, including QVs, into a MutationScorer, which is a
    # principled and fast way to test potential refinements to our
    # consensus sequence.
    r = cc.SparseSseQvRecursor()
    mms = cc.SparseSseQvMultiReadMutationScorer(r, PARAMS, roughConsensus)
    for cr in clippedReads[:QUIVER_COVERAGE]:
        mms.AddRead(AllQVsModel.extractFeatures(cr), int(cr.RCRefStrand))

    # To achieve greater accuracy (QV ~50, with this coverage level),
    # we iteratively refine the consensus sequence by identifying
    # mutations to the consensus that improve the likelihood (the
    # probability of the observed reads given the consensus as the
    # underlying template).
    css, didConverge = refineConsensus(mms)

    ga = cc.Align(refString[start:stop], mms.Template())
    print "  Quiver: %0.3f  (%d errors)" % (ga.Accuracy(), ga.Errors())
    table += Counter(ga.Transcript())
    print table
