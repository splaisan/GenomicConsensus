#################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

# Author: David Alexander

#
# This module makes the options globally available to all processes.
#
# Presently it relies on the fact that multiprocessing on Linux/Unix uses fork(),
# so the child processes inherit state from the main process.  If we ever wanted
# to port to Windows, we could support this module using the Manager protocol in
# in multiprocessing.
#
# Usage:
#  In the main process, before forking:
#    > from options import parseOptions, options
#    > parseOptions()
#    ...
#  then in any subprocess you can say
#    > from options import options
#  and get the loaded options dictionary.
#
from __future__ import absolute_import
import argparse, h5py, os, os.path, sys
from .utils import fileFormat
from . import __VERSION__

options = argparse.Namespace()

def consensusCoreVersion():
    try:
        import ConsensusCore
        return ConsensusCore.Version.VersionString()
    except:
        return None

def parseOptions(relax=False):
    """Parse the options and perform some due diligence on them, allowing for
    unit tests to relax things a bit.
    """
    desc = "Compute genomic consensus and call variants relative to the reference.."
    parser = argparse.ArgumentParser(description=desc)

    def checkInputFile(path):
        if not os.path.isfile(path):
            parser.error("Input file %s not found." % (path,))

    def checkOutputFile(path):
        try:
            f = open(path, "a")
            f.close()
        except:
            parser.error("Output file %s cannot be written." % (path,))

    parser.add_argument("inputFilename",
                        type=str,
                        help="The input cmp.h5 file")
    parser.add_argument("-o", "--outputFilename",
                        dest="outputFilenames",
                        required=not relax,
                        type=str,
                        action="append",
                        default=[],
                        help="The output filename(s), as a comma-separated list.")
    parser.add_argument("--verbose",
                        "-v",
                        dest="verbosity",
                        action="count",
                        help="Set the verbosity level.")
    parser.add_argument("-j",
                        "--numWorkers",
                        dest="numWorkers",
                        type=int,
                        default=1,
                        help="The number of worker processes to be used")
    parser.add_argument("--profile",
                        action="store_true",
                        dest="doProfiling",
                        default=False,
                        help="Enable Python-level profiling (using cProfile).")
    parser.add_argument("--referenceFilename", "--reference", "-r",
                        action="store",
                        dest="referenceFilename",
                        type=str,
                        required=not relax,
                        help="The filename of the reference FASTA file")

    # Since the reference isn't loaded at options processing time, we
    # can't grok the referenceWindow specified until later.  We store
    # it as a string (referenceWindowAsString) and it will later be
    # interpreted and stored as a proper window tuple (referenceWindow)
    parser.add_argument("--referenceWindow", "-w",
                        action="store",
                        dest="referenceWindowAsString",
                        type=str,
                        help="The window of the reference to be processed, in the format" + \
                             " refGroup:refStart-refEnd (default: entire reference).    "     )

    parser.add_argument("--algorithm",
                        action="store",
                        dest="algorithm",
                        type=str,
                        default="plurality")
    parser.add_argument("--parametersFile", "-P",
                        dest="parametersFile",
                        type=str,
                        default=None,
                        help="Parameter set filename (QuiverParameters.ini), or directory D " + \
                             "such that either D/*/GenomicConsensus/QuiverParameters.ini, "   + \
                             "or D/GenomicConsensus/QuiverParameters.ini, is found.  In the " + \
                             "former case, the lexically largest path is chosen.")
    parser.add_argument("--parameterSet", "-p",
                        action="store",
                        dest="parameterSet",
                        type=str,
                        default=None,
                        help="Name of parameter set to select from the parameters file.")

    parser.add_argument("--coverage", "-X",
                        action="store",
                        dest="coverage",
                        type=int,
                        default=100,
                        help="A designation of the maximum coverage level to be used for analysis." + \
                             " Exact interpretation is algorithm-specific.")
    parser.add_argument("--minMapQV", "-m",
                        action="store",
                        dest="minMapQV",
                        type=float,
                        default=10,
                        help="The minimum MapQV for reads that will be used for analysis.")
    parser.add_argument("--referenceChunkSize", "-C",
                        action="store",
                        dest="referenceChunkSize",
                        type=int,
                        default=500)
    parser.add_argument("--queueSize", "-Q",
                        action="store",
                        dest="queueSize",
                        type=int,
                        default=200)
    parser.add_argument("--threaded", "-T",
                        action="store_true",
                        dest="threaded",
                        default=False,
                        help="Run threads instead of processes (for debugging purposes only)")
    parser.add_argument("--dumpEvidence", "-d",
                        dest="dumpEvidence",
                        nargs="?",
                        default=None,
                        const="variants",
                        choices=["variants", "all"])
    parser.add_argument("--evidenceDirectory",
                        default="evidence_dump")

    def parseReadStratum(s):
        rs = map(int, s.split("/"))
        assert len(rs) == 2
        assert rs[0] < rs[1]
        return rs

    parser.add_argument(
        "--readStratum",
        help="A string of the form 'n/N', where n, and N are integers, 0 <= n < N, designating" \
             " that the reads are to be deterministically split into N strata of roughly even"  \
             " size, and stratum n is to be used for variant and consensus calling.  This is"   \
             " mostly useful for Quiver development.",
        dest="readStratum",
        default=None,
        type=parseReadStratum)

    parser.add_argument("--minConfidence", "-q",
                        action="store",
                        dest="minConfidence",
                        type=int,
                        default=40,
                        help="The minimum confidence for a variant call to be output to variants.gff")
    parser.add_argument("--minCoverage", "-x",
                        action="store",
                        dest="minCoverage",
                        default=5,
                        type=int,
                        help="The minimum site coverage that must be achieved for variant calls and " + \
                             "consensus to be calculated for a site.")
    parser.add_argument("--referenceChunkOverlap",
                        action="store",
                        dest="referenceChunkOverlap",
                        type=int)
    parser.add_argument(
        "--noEvidenceConsensusCall",
        action="store",
        choices=["nocall", "reference", "lowercasereference"],
        default="lowercasereference",
        help="The consensus base that will be output for sites with no effective coverage.")
    parser.add_argument(
        "--aligner", "-a",
        action="store",
        choices=["affine", "simple"],
        default="affine",
        help="The pairwise alignment algorithm that will be used to produce variant calls" \
             " from the consensus (Quiver only).")
    parser.add_argument(
        "--refineDinucleotideRepeats",
        default=True,
        help="Require quiver maximum likelihood search to try one less/more repeat copy in"  \
             " dinucleotide repeats, which seem to be the most frequent cause of suboptimal" \
             " convergence (getting trapped in local optimum) (Quiver only)")

    parser.add_argument(
        "--fancyChunking",
        default=False,
        action="store_true",
        help="Currently experimental chunking designed to handle coverage cutouts better")

    class PrintVersionAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print "  GenomicConsensus version: %s" % __VERSION__
            print "  ConsensusCore version: %s" % \
                (consensusCoreVersion() or "ConsensusCore unavailable")
            print "  h5py version: %s" % h5py.version.version
            print "  hdf5 version: %s" % h5py.version.hdf5_version
            sys.exit(0)

    parser.add_argument("--version",
                        nargs=0,
                        action=PrintVersionAction)

    parser.parse_args(namespace=options)

    options.gffOutputFilename   = None
    options.fastaOutputFilename = None
    options.fastqOutputFilename = None
    options.csvOutputFilename   = None

    for outputFilename in options.outputFilenames:
        fmt = fileFormat(outputFilename)
        if   fmt == "GFF":   options.gffOutputFilename   = outputFilename
        elif fmt == "FASTA": options.fastaOutputFilename = outputFilename
        elif fmt == "FASTQ": options.fastqOutputFilename = outputFilename
        elif fmt == "CSV":   options.csvOutputFilename   = outputFilename

    # chill-out, stop worrying so much
    if relax: return

    for path in (options.inputFilename, options.referenceFilename):
        if path != None:
            checkInputFile(path)

    for path in options.outputFilenames:
        if path != None:
            checkOutputFile(path)

def importAdditionalDefaultOptions(additionalDefaults):
    # After parsing the arguments, we need to patch the options to
    # include algorithm-specific defaults, if they have not been
    # overriden on the command-line.
    optionsDictView = vars(options)
    for k, v in additionalDefaults.iteritems():
        if k not in optionsDictView or optionsDictView[k] is None:
            optionsDictView[k] = v
