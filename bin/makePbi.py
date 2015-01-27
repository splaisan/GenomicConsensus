#!/usr/bin/env python

import argparse
import h5py
import pysam
import numpy as np
from collections import Counter, OrderedDict

from pbcore.io import BamReader, BamAlignment

# Call: makePbi.py bamfile referenceFasta
# Warning: this is a very inefficient reference implementation.

PBI_VERSION = "3.0"

PBI_COLUMNS_AND_TYPES = [ ("tId"               , np.int32),
                          ("tStart"            , np.int32),
                          ("tEnd"              , np.int32),
                          ("qId"               , np.int32),
                          ("qStart"            , np.int32),
                          ("qEnd"              , np.int32),
                          ("rStart"            , np.int32),
                          ("rEnd"              , np.int32),
                          ("holeNumber"        , np.int32),
                          ("isReverseStrand"   , np.int8),
                          ("nM"                , np.int32),
                          ("nMM"               , np.int32),
                          ("mapQV"             , np.uint8),
                          ("virtualFileOffset" , np.uint64) ]

def main():
    parser = argparse.ArgumentParser(description="Build PacBio BAM index for a bam")
    parser.add_argument("bamFile", type=argparse.FileType("r"))
    parser.add_argument("--referenceFasta", type=argparse.FileType("r"))
    parser.add_argument("--lite", action="store_true")

    args = parser.parse_args()

    bamFname = args.bamFile.name
    if args.referenceFasta:
        refFname = args.referenceFasta.name
    else:
        refFname = None

    pbi = h5py.File(bamFname + ".pbi", "w")
    p = pysam.Samfile(bamFname, check_sq=False)
    B = BamReader(bamFname, refFname)

    lsts = dict((columnName, [])
                for (columnName, dtype_) in PBI_COLUMNS_AND_TYPES)
    p.reset()
    while True:
        offset = p.tell()
        try:
            rawAln = next(p)
        except StopIteration:
            break

        aln = BamAlignment(B, rawAln)

        if aln.isMapped:
            lsts["tId"               ].append(aln.tId)
            lsts["tStart"            ].append(aln.tStart)
            lsts["tEnd"              ].append(aln.tEnd)
            lsts["qId"               ].append(aln.qId)
            lsts["qStart"            ].append(aln.qStart)
            lsts["qEnd"              ].append(aln.qEnd)
            lsts["rStart"            ].append(aln.rStart)
            lsts["rEnd"              ].append(aln.rEnd)
            lsts["holeNumber"        ].append(aln.HoleNumber)
            lsts["isReverseStrand"   ].append(aln.isReverseStrand)
            lsts["mapQV"             ].append(aln.MapQV)
            lsts["virtualFileOffset" ].append(offset)
            if not args.lite:
                transcript = aln.transcript()
                moveCounts = Counter(transcript)
                lsts["nM"            ].append(moveCounts["M"])
                lsts["nMM"           ].append(moveCounts["R"])

        else:
            # Unmapped
            lsts["tId"               ].append(-1)
            lsts["tStart"            ].append(-1)
            lsts["tEnd"              ].append(-1)
            lsts["qId"               ].append(aln.qId)
            lsts["qStart"            ].append(aln.qStart)
            lsts["qEnd"              ].append(aln.qEnd)
            lsts["rStart"            ].append(-1)
            lsts["rEnd"              ].append(-1)
            lsts["holeNumber"        ].append(aln.HoleNumber)
            lsts["isReverseStrand"   ].append(-1)
            lsts["mapQV"             ].append(-1)
            lsts["virtualFileOffset" ].append(offset)
            lsts["nM"                ].append(-1)
            lsts["nMM"               ].append(-1)


    # Lists to arrays
    dsets = {}
    for (name, dtype) in PBI_COLUMNS_AND_TYPES:
        dsets[name] = np.fromiter(lsts[name], dtype)
        del lsts[name]

    # Write to file, gzipped
    grp = pbi.create_group("PacBioBamIndex")
    grp.attrs["Version"] = np.string_(PBI_VERSION)
    columnsGrp = grp.create_group("Columns")
    for (name, ds) in dsets.iteritems():
        columnsGrp.create_dataset(name, data=ds, chunks=True, compression="gzip")
    pbi.close()



if __name__ == '__main__':
    main()
