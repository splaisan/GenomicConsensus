#!/usr/bin/env python

#
# A tool for viewing the alignments in a given reference window.
# Will be migrated to pbh5tools when possible.
#

import argparse, sys, numpy as np
from pbcore.io import CmpH5Reader

def formatReferenceCoordinates(refId, refStart, refEnd):
    canvas = np.zeros(refEnd - refStart, dtype="S1")
    canvas[:] = " "
    refCoords = np.arange(refStart, refEnd)
    for refPos in refCoords:
        if refPos % 20 == 0:
            sPos = str(refPos)
            sPosLen = len(sPos)
            cursorPos = refPos - sPosLen - refStart + 1
            if cursorPos > 0:
                canvas[cursorPos:(cursorPos+sPosLen)] = \
                    np.fromstring(sPos, dtype="S1")
    return canvas.tostring()

def formatSeparatorLine(refId, refStart, refEnd):
    # Formats line like "--|----+----|----+--"
    canvas = np.zeros(refEnd - refStart, dtype="S1")
    canvasCoords = np.arange(refStart, refEnd)
    canvas[:] = "="
    canvas[canvasCoords % 10 == 0] = "|"
    canvas[canvasCoords % 10 == 5] = "+"
    return canvas.tostring()

def formatReference(refWindow, referenceByIdDict):
    pass

def formatAlignedRead(refId, refStart, refEnd, cmpH5, rowNumber):
    assert refEnd > refStart
    canvas = np.zeros(refEnd - refStart, dtype="S1")
    canvas[:] = " "
    clippedRead = cmpH5[rowNumber].clippedTo(refStart, refEnd)
    fullRead = np.fromstring(clippedRead.read(orientation="genomic"), dtype="S1")
    transcript = np.fromstring(clippedRead.transcript(orientation="genomic"), dtype="S1")
    noInsertionsRead = np.extract(transcript != "I", fullRead)
    canvas[(clippedRead.tStart-refStart):(clippedRead.tEnd-refStart)] = noInsertionsRead
    return canvas.tostring()

def formatAlignedReads(refWindow, cmpH5, rowNumbers):
    refId, refStart, refEnd = refWindow
    return [ formatAlignedRead(refId, refStart, refEnd, cmpH5, rowNumber)
             for rowNumber in rowNumbers ]

def main():
    parser = argparse.ArgumentParser(description="View alignments")
    parser.add_argument("cmpH5Filename", type=str)
    parser.add_argument("refId", type=str)
    parser.add_argument("refStart", type=int)
    parser.add_argument("refEnd", type=int)
    parser.add_argument("--depth", "-X", type=int, default=20)
    parser.add_argument("--minMapQV", type=int, default=10)
    parser.add_argument("--rowNumbers", type=int, nargs="+", default=None)
    parser.add_argument("--columns", type=str, nargs="+", default=None)

    options = parser.parse_args()

    cmpH5 = CmpH5Reader(options.cmpH5Filename)
    refId = cmpH5.referenceInfo(options.refId).ID
    refWindow = (refId, options.refStart, options.refEnd)

    if options.rowNumbers != None:
        rowNumbers = options.rowNumbers
    else:
        rowNumbers = [ rowNumber
                       for rowNumber in cmpH5.readsInRange(*refWindow, justIndices=True)
                       if cmpH5.MapQV[rowNumber] > options.minMapQV
                     ][:options.depth]

    preMargin = " " * 10
    print preMargin + formatReferenceCoordinates(*refWindow)
    print preMargin + formatSeparatorLine(*refWindow)
    for rn, ar in zip(rowNumbers, formatAlignedReads(refWindow, cmpH5, rowNumbers)):
        print ("%8d  " % rn)  + ar


if __name__ == '__main__':
    main()
