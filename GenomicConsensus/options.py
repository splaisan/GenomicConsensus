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
import argparse, os, os.path, sys
from .utils import fileFormat

options = argparse.Namespace()

VERSION = "0.2.0"

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
    parser.add_argument("--parameters", "-p",
                        action="store",
                        dest="parameters",
                        type=str,
                        default=None)

    parser.add_argument("--coverage", "-X",
                        action="store",
                        dest="coverage",
                        type=int,
                        default=None,
                        help="A designation of the maximum coverage level to be used for analysis." + \
                             " Exact interpretation is algorithm-specific.")
    parser.add_argument("--mapQvThreshold", "-m",
                        action="store",
                        dest="mapQvThreshold",
                        type=float,
                        default=10.0)
    parser.add_argument("--referenceChunkSize", "-C",
                        action="store",
                        dest="referenceChunkSize",
                        type=int,
                        default=1000)
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
    #
    # Default values for these arguments are algorithm specific.
    # See 'additionalDefaultOptions' in each algorithm module.
    #
    parser.add_argument("--variantConfidenceThreshold", "-q",
                        action="store",
                        dest="variantConfidenceThreshold",
                        type=float)
    parser.add_argument("--variantCoverageThreshold", "-x",
                        action="store",
                        dest="variantCoverageThreshold",
                        type=int)
    parser.add_argument("--referenceChunkOverlap",
                        action="store",
                        dest="referenceChunkOverlap",
                        type=int)
    parser.add_argument("--noEvidenceConsensusCall",
                        action="store",
                        choices=["nocall", "reference"],
                        default="nocall",
                        help="The consensus base that will be output for sites with no effective coverage.")
    parser.add_argument("--aligner", "-a",
                        action="store",
                        choices=["affine", "simple"],
                        default="affine",
                        help="The pairwise alignment algorithm that will be used to produce variant calls" \
                             " from the consensus (Quiver only).")

    class PrintVersionAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            print "  Genomic Consensus version %s" % VERSION
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

def importAdditionalDefaulOptions(additionalDefaults):
    # After parsing the arguments, we need to patch the options to
    # include algorithm-specific defaults, if they have not been
    # overriden on the command-line.
    optionsDictView = vars(options)
    for k, v in additionalDefaults.iteritems():
        if k not in optionsDictView or optionsDictView[k] == None:
            optionsDictView[k] = v
