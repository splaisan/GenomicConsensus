#!/usr/bin/env python

import argparse, gzip, numpy as np, sys
from collections import namedtuple

from pbcore.io import GffReader, GffWriter, Gff3Record, parseGffLine

consensusCsvDtypes =  [("referenceId"         , np.uint32),
                       ("referencePos"        , np.uint32),
                       ("coverage"            , np.uint32),
                       ("consensus"           , np.object),
                       ("consensusConfidence" , np.uint8),
                       ("consensusFrequency"  , np.uint32)]

#
# Note: GFF-style coordinates
#
Region = namedtuple("Region", ("seqid", "start", "end"))


def lookup(key, alist):
    for k, v in alist:
        if k == key:
            return v

def refNameToId(s):
    # remove 'ref'
    return int(s[3:])

def main():
    headers = [
        ("source", "GenomicConsensus 0.2.0"),
        ("pacbio-alignment-summary-version", "0.5"),
        ("source-commandline", " ".join(sys.argv)),
        ]

    desc = "Augment the alignment_summary.gff file with consensus and variants information."
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument("--variantsGff",
                        type=str,
                        help="Input variants.gff or variants.gff.gz filename",
                        required=True)
    parser.add_argument("--consensusCsv",
                        type=str,
                        help="Input consensus.csv or consensus.csv.gz filename",
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
    if options.consensusCsv.endswith(".gz"):
        f = gzip.open(options.consensusCsv)
    else:
        f = open(options.consensusCsv)
    inputConsensusCsv = np.loadtxt(f, dtype=consensusCsvDtypes, delimiter=",",
                                   skiprows=1, converters={0 : refNameToId})
    inputConsensusCsv = inputConsensusCsv.view(np.recarray)


    inputAlignmentSummaryGff = GffReader(options.inputAlignmentSummaryGff)

    summaries = {}
    for gffRecord in inputAlignmentSummaryGff:
        region = Region(refNameToId(gffRecord.seqid), gffRecord.start, gffRecord.end)
        summaries[region] = { "ins" : 0,
                              "del" : 0,
                              "sub" : 0,
                              "cQv" : (0, 0, 0)
                             }
    inputAlignmentSummaryGff.close()

    for region in summaries:
        regionTbl = inputConsensusCsv[(inputConsensusCsv.referenceId == region.seqid) &
                                      (region.start-1 <= inputConsensusCsv.referencePos) &
                                      (inputConsensusCsv.referencePos < region.end)]
        regionQuality = regionTbl.consensusConfidence
        if len(regionQuality) > 0:
            summaries[region]["cQv"] = (np.min(regionQuality),
                                        np.median(regionQuality),
                                        np.max(regionQuality))

    counterNames = { "insertion"    : "ins",
                     "deletion"     : "del",
                     "substitution" : "sub" }
    for variantGffRecord in inputVariantsGff:
        variantGffRecordRefId = refNameToId(variantGffRecord.seqid)
        for region in summaries:
            if (region.seqid == variantGffRecordRefId and
                region.start <= variantGffRecord.start <= region.end):
                summary = summaries[region]
                counterName = counterNames[variantGffRecord.type]
                variantLength = int(lookup("length", variantGffRecord.attributes))
                summary[counterName] += variantLength

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
            summary = summaries[Region(refNameToId(rec.seqid), rec.start, rec.end)]
            if "cQv" in summary:
                cQvTuple = summary["cQv"]
                line += ";%s=%s" % ("cQv", ",".join(str(int(f)) for f in cQvTuple))
            for counterName in counterNames.values():
                if counterName in summary:
                    line += ";%s=%d" % (counterName, summary[counterName])
            print >>outputAlignmentSummaryGff, line

if __name__ == "__main__":
    main()


