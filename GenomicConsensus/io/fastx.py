from __future__ import absolute_import
import os, numpy as np

OUTPUT_CHUNK_SIZE = 10**6

def qvArrayToString(qvArray):
    assert qvArray.dtype == np.uint8
    return np.minimum(33 + qvArray, 126).tostring()

class FastqWriter(object):
    """
    For encoding details see:
    http://en.wikipedia.org/wiki/FASTQ_format#Encoding
    """
    def __init__(self, file):
        self._file = file

    def writeRecord(self, sequenceName, sequenceArray, qvArray):
        self._file.write("@" + sequenceName + os.linesep)

        # write a chunk at a time, so we don't try to allocate
        # massive temporary strings / arrays.
        for offset in xrange(0, len(sequenceArray), OUTPUT_CHUNK_SIZE):
            self._file.write("".join(sequenceArray[offset:(offset + OUTPUT_CHUNK_SIZE)]))

        self._file.write(os.linesep)
        self._file.write("+" + os.linesep)

        for offset in xrange(0, len(sequenceArray), OUTPUT_CHUNK_SIZE):
            # we need the qv values in consensus coordinates.
            consensusLens = [len(s)
                             for s in sequenceArray[offset:(offset + OUTPUT_CHUNK_SIZE)]]
            consensusQvs = np.repeat(qvArray[offset:(offset + OUTPUT_CHUNK_SIZE)],
                                     consensusLens)
            self._file.write(qvArrayToString(consensusQvs))

        self._file.write(os.linesep)


class FastaWriter(object):
    def __init__(self, file):
        self._file = file

    def writeRecord(self, sequenceName, sequenceArray):
        self._file.write(">" + sequenceName + os.linesep)
        for offset in xrange(0, len(sequenceArray), OUTPUT_CHUNK_SIZE):
            self._file.write("".join(sequenceArray[offset:(offset + OUTPUT_CHUNK_SIZE)]))
        self._file.write(os.linesep)
