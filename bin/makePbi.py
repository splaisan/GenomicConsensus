#!/usr/bin/env python

import argparse
import h5py
import pysam
import numpy as np
from collections import Counter, OrderedDict

from GenomicConsensus.io import BamReader, BamAlignment

# Call: makePbi.py bamfile referenceFasta

PBI_VERSION = "0.1"

PBI_COLUMNS_AND_TYPES = [ ("ReadGroupID"       , np.uint32),
                          ("tId"               , np.uint32),
                          ("tStart"            , np.uint32),
                          ("tEnd"              , np.uint32),
                          ("isReverseStrand"   , np.uint8),
                          ("HoleNumber"        , np.uint32),
                          ("rStart"            , np.uint32),
                          ("rEnd"              , np.uint32),
                          ("MapQV"             , np.uint8),
                          ("nM"                , np.uint32),
                          ("nMM"               , np.uint32),
                          ("nIns"              , np.uint32),
                          ("nDel"              , np.uint32),
                          ("virtualFileOffset" , np.uint64)  ]

def main():
    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument("bamFile", type=argparse.FileType("r"))
    parser.add_argument("referenceFasta", type=argparse.FileType("r"))
    args = parser.parse_args()

    bamFname = args.bamFile.name
    refFname = args.referenceFasta.name

    pbi = h5py.File(bamFname + ".pbi", "w")

    B = BamReader(bamFname, refFname)
    p = pysam.Samfile(bamFname)

    shape = (len(B),)
    #shape = (len(B), len(PBI_COLUMNS_AND_TYPES))

    dsets = OrderedDict((columnName, np.zeros(shape=shape, dtype=dtype_))
                        for (columnName, dtype_) in PBI_COLUMNS_AND_TYPES)
    p.reset()
    for i in xrange(len(B)):
        offset = p.tell()
        rawAln = next(p)
        aln = BamAlignment(B, rawAln)

        transcript = aln.transcript()
        moveCounts = Counter(transcript)

        dsets["ReadGroupID"]      [i] = aln.readGroup.ID
        dsets["tId"]              [i] = aln.tId
        dsets["tStart"]           [i] = aln.tStart
        dsets["tEnd"]             [i] = aln.tEnd
        dsets["isReverseStrand"]  [i] = aln.isReverseStrand
        dsets["HoleNumber"]       [i] = aln.HoleNumber
        dsets["rStart"]           [i] = aln.rStart
        dsets["rEnd"]             [i] = aln.rEnd
        dsets["MapQV"]            [i] = aln.MapQV
        dsets["nM"]               [i] = moveCounts["M"]
        dsets["nMM"]              [i] = moveCounts["R"]
        dsets["nIns"]             [i] = moveCounts["I"]
        dsets["nDel"]             [i] = moveCounts["D"]
        dsets["virtualFileOffset"][i] = offset

    # Write to file, gzipped

    grp = pbi.create_group("PacBioBamIndex")
    grp.attrs["Version"] = np.string_(PBI_VERSION)
    columnsGrp = grp.create_group("Columns")
    for (name, ds) in dsets.iteritems():
        columnsGrp.create_dataset(name, data=ds, chunks=True, compression="gzip")
    pbi.close()



if __name__ == '__main__':
    main()
