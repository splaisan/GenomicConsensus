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

def parseOptions():
    """
    Parse the options and perform some due diligence on them
    """
    desc = "Compute genomic consensus and call variants relative to the reference."
    parser = argparse.ArgumentParser(description=desc, add_help=False)

    def canonicalizedFilePath(path):
        return os.path.abspath(os.path.expanduser(path))

    def checkInputFile(path):
        if not os.path.isfile(path):
            parser.error("Input file %s not found." % (path,))

    def checkOutputFile(path):
        try:
            f = open(path, "a")
            f.close()
        except:
            parser.error("Output file %s cannot be written." % (path,))


    basics = parser.add_argument_group("Basic required options")
    basics.add_argument(
        "inputFilename",
        type=canonicalizedFilePath,
        help="The input cmp.h5 file")
    basics.add_argument(
        "--referenceFilename", "--reference", "-r",
        action="store",
        dest="referenceFilename",
        type=canonicalizedFilePath,
        required=True,
        help="The filename of the reference FASTA file")
    basics.add_argument(
        "-o", "--outputFilename",
        dest="outputFilenames",
        required=True,
        type=str,
        action="append",
        default=[],
        help="The output filename(s), as a comma-separated list." + \
             "Valid output formats are .fa/.fasta, .fq/.fastq, .gff")

    parallelism = parser.add_argument_group("Parallelism")
    parallelism.add_argument(
        "-j", "--numWorkers",
        dest="numWorkers",
        type=int,
        default=1,
        help="The number of worker processes to be used")

    filtering = parser.add_argument_group("Output filtering")
    filtering.add_argument(
        "--minConfidence", "-q",
        action="store",
        dest="minConfidence",
        type=int,
        default=40,
        help="The minimum confidence for a variant call to be output to variants.gff")
    filtering.add_argument(
        "--minCoverage", "-x",
        action="store",
        dest="minCoverage",
        default=5,
        type=int,
        help="The minimum site coverage that must be achieved for variant calls and " + \
             "consensus to be calculated for a site.")
    filtering.add_argument(
        "--noEvidenceConsensusCall",
        action="store",
        choices=["nocall", "reference", "lowercasereference"],
        default="lowercasereference",
        help="The consensus base that will be output for sites with no effective coverage.")


    readSelection = parser.add_argument_group("Read selection/filtering")
    readSelection.add_argument(
        "--coverage", "-X",
        action="store",
        dest="coverage",
        type=int,
        default=100,
        help="A designation of the maximum coverage level to be used for analysis." + \
             " Exact interpretation is algorithm-specific.")
    readSelection.add_argument(
        "--minMapQV", "-m",
        action="store",
        dest="minMapQV",
        type=float,
        default=10,
        help="The minimum MapQV for reads that will be used for analysis.")
    # Since the reference isn't loaded at options processing time, we
    # can't grok the referenceWindow specified until later.  We store
    # it as a string (referenceWindowsAsString) and it will later be
    # interpreted and stored as a proper window tuple (referenceWindow)
    readSelection.add_argument(
        "--referenceWindow", "--referenceWindows", "-w",
        action="store",
        dest="referenceWindowsAsString",
        type=str,
        help="The window (or multiple comma-delimited windows) of the reference to " + \
             "be processed, in the format refGroup:refStart-refEnd "                 + \
             "(default: entire reference).",
        default=None)

    def slurpWindowFile(fname):
        return ",".join(map(str.strip, open(fname).readlines()))

    readSelection.add_argument(
        "--referenceWindowsFile", "-W",
        action="store",
        dest="referenceWindowsAsString",
        type=slurpWindowFile,
        help="A file containing reference window designations, one per line",
        default=None)
    readSelection.add_argument(
        "--barcode",
        type=str,
        dest="_barcode",
        help="Only process reads with the given barcode name.")
    def parseReadStratum(s):
        rs = map(int, s.split("/"))
        assert len(rs) == 2
        assert rs[0] < rs[1]
        return rs
    readSelection.add_argument(
        "--readStratum",
        help="A string of the form 'n/N', where n, and N are integers, 0 <= n < N, designating" \
             " that the reads are to be deterministically split into N strata of roughly even"  \
             " size, and stratum n is to be used for variant and consensus calling.  This is"   \
             " mostly useful for Quiver development.",
        dest="readStratum",
        default=None,
        type=parseReadStratum)


    algorithm = parser.add_argument_group("Algorithm and parameter settings")
    algorithm.add_argument(
        "--algorithm",
        action="store",
        dest="algorithm",
        type=str,
        default="quiver")
    algorithm.add_argument(
        "--parametersFile", "-P",
        dest="parametersFile",
        type=str,
        default=None,
        help="Parameter set filename (QuiverParameters.ini), or directory D " + \
             "such that either D/*/GenomicConsensus/QuiverParameters.ini, "   + \
             "or D/GenomicConsensus/QuiverParameters.ini, is found.  In the " + \
             "former case, the lexically largest path is chosen.")
    algorithm.add_argument(
        "--parametersSpec", "-p",
        action="store",
        dest="parametersSpec",
        type=str,
        default="auto",
        help="Name of parameter set (chemistry.model) to select from the "   + \
             "parameters file, or just the name of the chemistry, in which " + \
             "case the best available model is chosen.  Default is 'auto', " + \
             "which selects the best parameter set from the cmp.h5")

    debugging = parser.add_argument_group("Verbosity and debugging/profiling")
    debugging.add_argument("--help", "-h",
                           action="help")
    class PrintVersionAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print "  GenomicConsensus version: %s" % __VERSION__
            print "  ConsensusCore version: %s" % \
                (consensusCoreVersion() or "ConsensusCore unavailable")
            print "  h5py version: %s" % h5py.version.version
            print "  hdf5 version: %s" % h5py.version.hdf5_version
            sys.exit(0)
    debugging.add_argument("--version",
                           nargs=0,
                           action=PrintVersionAction)
    debugging.add_argument(
        "--verbose", "-v",
        dest="verbosity",
        action="count",
        help="Set the verbosity level.")
    debugging.add_argument(
        "--quiet",
        dest="quiet",
        action="store_true",
        help="Turn off all logging, including warnings")
    debugging.add_argument(
        "--debug",
        action="store_true",
        dest="doDebugging",
        default=False,
        help="Enable Python-level debugging (using pdb).")
    debugging.add_argument(
        "--profile",
        action="store_true",
        dest="doProfiling",
        default=False,
        help="Enable Python-level profiling (using cProfile).")
    debugging.add_argument(
        "--dumpEvidence", "-d",
        dest="dumpEvidence",
        nargs="?",
        default=None,
        const="variants",
        choices=["variants", "all"])
    debugging.add_argument(
        "--evidenceDirectory",
        default="evidence_dump")
    debugging.add_argument(
        "--annotateGFF",
        action="store_true",
        help="Augment GFF variant records with additional information")

    advanced = parser.add_argument_group("Advanced configuration options")
    advanced.add_argument(
        "--diploid",
        action="store_true",
        help="Enable detection of heterozygous variants (experimental)")
    advanced.add_argument(
        "--queueSize", "-Q",
        action="store",
        dest="queueSize",
        type=int,
        default=200)
    advanced.add_argument(
        "--threaded", "-T",
        action="store_true",
        dest="threaded",
        default=False,
        help="Run threads instead of processes (for debugging purposes only)")
    advanced.add_argument(
        "--referenceChunkSize", "-C",
        action="store",
        dest="referenceChunkSize",
        type=int,
        default=500)
    advanced.add_argument(
        "--fancyChunking",
        default=True,
        action="store_true",
        help="Adaptive reference chunking designed to handle coverage cutouts better")
    advanced.add_argument(
        "--simpleChunking",
        dest="fancyChunking",
        action="store_false",
        help="Disable adaptive reference chunking")
    advanced.add_argument(
        "--referenceChunkOverlap",
        action="store",
        dest="referenceChunkOverlap",
        type=int,
        default=5)
    advanced.add_argument(
        "--autoDisableHdf5ChunkCache",
        action="store",
        type=int,
        default=500,
        help="Disable the HDF5 chunk cache when the number of datasets in the cmp.h5 " + \
             "exceeds the given threshold")
    advanced.add_argument(
        "--aligner", "-a",
        action="store",
        choices=["affine", "simple"],
        default="affine",
        help="The pairwise alignment algorithm that will be used to produce variant calls" \
             " from the consensus (Quiver only).")
    advanced.add_argument(
        "--refineDinucleotideRepeats",
        dest="refineDinucleotideRepeats",
        action="store_true",
        help="Require quiver maximum likelihood search to try one less/more repeat copy in"  \
             " dinucleotide repeats, which seem to be the most frequent cause of suboptimal" \
             " convergence (getting trapped in local optimum) (Quiver only)")
    advanced.add_argument(
        "--noRefineDinucleotideRepeats",
        dest="refineDinucleotideRepeats",
        action="store_false",
        help="Disable dinucleotide refinement")
    advanced.set_defaults(refineDinucleotideRepeats=True)
    advanced.add_argument(
        "--fast",
        dest="fastMode",
        action="store_true",
        help="Cut some corners to run faster.  Unsupported!")
    advanced.add_argument(
        "--skipUnrecognizedContigs",
        action="store_true",
        help="Do not abort when told to process a reference window (via -w/--referenceWindow[s]) " \
             "that has no aligned coverage.  Outputs emptyish files if there are no remaining "    \
             "non-degenerate windows.  Only intended for use by smrtpipe scatter/gather.")

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

    for path in (options.inputFilename, options.referenceFilename):
        if path != None:
            checkInputFile(path)

    for path in options.outputFilenames:
        if path != None:
            checkOutputFile(path)

    options.shellCommand = " ".join(sys.argv)


def resolveOptions(cmpH5):
    """
    Some of the options are provided as strings by the user, but need
    to be translated into internal identifiers.  These options are
    encoded as options._optionName; here we lookup the ID and store it
    as options.optionName.

    This is essentially just an order-of-initialization issue.
    """
    if options._barcode != None:
        if not cmpH5.isBarcoded:
            raise Exception("cmp.h5 file is not barcoded!")
        if options._barcode not in cmpH5.barcode:
            raise Exception("Barcode with given name not present in cmp.h5 file!")
        options.barcode = cmpH5.barcode[options._barcode]
    else:
        options.barcode = None
