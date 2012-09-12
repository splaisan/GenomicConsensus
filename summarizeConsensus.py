#!/usr/bin/env python

import argparse, gzip, numpy as np, sys
from collections import namedtuple

from pbcore.io import GffReader, GffWriter, Gff3Record, parseGffLine
from GenomicConsensus.utils import error_probability_to_qv

#
# Note: GFF-style coordinates
#
Region = namedtuple("Region", ("seqid", "start", "end"))


def lookup(key, alist):
    for k, v in alist:
        if k == key:
            return v

def main():
    headers = [
        ("source", "GenomicConsensus 0.2.0"),
        ("pacbio-alignment-summary-version", "0.6"),
        ("source-commandline", " ".join(sys.argv)),
        ]

    desc = "Augment the alignment_summary.gff file with consensus and variants information."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--variantsGff",
                        type=str,
                        help="Input variants.gff or variants.gff.gz filename",
                        required=True)
    parser.add_argument("--output",
                        "-o",
                        type=str,
                        help="Output alignment_summary.gff filename")
    parser.add_argument("inputAlignmentSummaryGff",
                        type=str,
                        help="Input alignment_summary.gff filename")

    options = parser.parse_args()

    inputVariantsGff = GffReader(options.variantsGff)
    inputAlignmentSummaryGff = GffReader(options.inputAlignmentSummaryGff)

    summaries = {}
    for gffRecord in inputAlignmentSummaryGff:
        region = Region(gffRecord.seqid, gffRecord.start, gffRecord.end)
        summaries[region] = { "ins" : 0,
                              "del" : 0,
                              "sub" : 0,
                              "cQv" : (0, 0, 0)
                             }
    inputAlignmentSummaryGff.close()

    counterNames = { "insertion"    : "ins",
                     "deletion"     : "del",
                     "substitution" : "sub" }
    for variantGffRecord in inputVariantsGff:
        for region in summaries:
            summary = summaries[region]
            if (region.seqid == variantGffRecord.seqid and
                region.start <= variantGffRecord.start <= region.end):
                counterName = counterNames[variantGffRecord.type]
                variantLength = int(lookup("length", variantGffRecord.attributes))
                summary[counterName] += variantLength
            # TODO: base consensusQV on effective coverage
            summary["cQv"] = (20, 20, 20)

    inputAlignmentSummaryGff = open(options.inputAlignmentSummaryGff)
    outputAlignmentSummaryGff = open(options.output, "w")

    inHeader = True

    for line in inputAlignmentSummaryGff:
        line = line.rstrip()

        # Pass any metadata line straight through
        if line[0] == "#":
            print >>outputAlignmentSummaryGff, line.strip()
            continue

        if inHeader:
            # We are at the end of the header -- write the tool-specific headers
            for k, v in headers:
                print >>outputAlignmentSummaryGff, ("##%s %s" % (k, v))
            inHeader = False

        # Parse the line
        rec = parseGffLine(line)

        if rec.type == "region":
            summary = summaries[(rec.seqid, rec.start, rec.end)]
            if "cQv" in summary:
                cQvTuple = summary["cQv"]
                line += ";%s=%s" % ("cQv", ",".join(str(int(f)) for f in cQvTuple))
            for counterName in counterNames.values():
                if counterName in summary:
                    line += ";%s=%d" % (counterName, summary[counterName])
            print >>outputAlignmentSummaryGff, line

if __name__ == "__main__":
    main()
