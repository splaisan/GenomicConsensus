from __future__ import absolute_import
import gzip, os, numpy as np, itertools

__all__ = [ "FastaWriter",
            "FastqWriter" ]

OUTPUT_CHUNK_SIZE = 10**6
LINE_WIDTH = 60

def chunks(n, l):
    for i in xrange(0, len(l), n):
        yield l[i:i+n]

def isFileLikeObject(o):
    return hasattr(o, "read") and hasattr(o, "write")

def getFileHandle(filenameOrFile, mode):
    """
    Given a filename not ending in ".gz", open the file with the
    appropriate mode.

    Given a filename ending in ".gz", return a filehandle to the
    unzipped stream.

    Given a file object, return it unless the mode is incorrect--in
    that case, raise an exception.
    """
    assert mode in ("r", "w")

    if isinstance(filenameOrFile, str):
        if filenameOrFile.endswith(".gz"):
            return gzip.open(filenameOrFile, mode)
        else:
            return open(filenameOrFile, mode)
    elif isFileLikeObject(filenameOrFile):
        return filenameOrFile
    else:
        raise Exception("Invalid type to getFileHandle")

class WrappedWriter(object):

    def __init__(self, fileObj, lineWidth):
        self.file = getFileHandle(fileObj, "w")
        self.lineWidth = lineWidth
        self._currentColumn = 0

    def newline(self):
        if self._currentColumn:
            self.file.write("\n")
            self._currentColumn = 0

    def write(self, contents):
        assert "\n" not in contents
        firstLineLen = self.lineWidth - self._currentColumn
        firstLine = contents[:firstLineLen]
        remainingLines = chunks(self.lineWidth, contents[firstLineLen:])
        self.file.write(firstLine)
        self._currentColumn += len(firstLine)
        for line in remainingLines:
            self.file.write("\n")
            self.file.write(line)
            self._currentColumn = len(line)

    def writeUnwrapped(self, contents):
        assert "\n" not in contents
        self.file.write(contents)
        self._currentColumn += len(contents)

    def close(self):
        if not self.file.closed:
            self.newline()
            self.file.close()

    def __del__(self):
        self.close()


class FastaWriter(object):

    def __init__(self, fileObj):
        self.wrappedWriter = WrappedWriter(fileObj, 60)

    def writeRecord(self, sequenceHeader, sequenceArray):
        """
        Sequence array is understood to be a string, or a numpy array
        containing strings which must be concatenated.
        """
        self.wrappedWriter.newline()
        self.wrappedWriter.writeUnwrapped(">" + sequenceHeader)
        self.wrappedWriter.newline()
        for chunk in chunks(OUTPUT_CHUNK_SIZE, sequenceArray):
            self.wrappedWriter.write("".join(chunk))


    def close(self):
        self.wrappedWriter.close()

    def __del__(self):
        self.close()

class FastqWriter(object):
    """
    For encoding details see:
    http://en.wikipedia.org/wiki/FASTQ_format#Encoding
    """
    def __init__(self, fileObj):
        self.file = getFileHandle(fileObj, "w")

    @staticmethod
    def _qvArrayToString(qvArray):
        assert qvArray.dtype == np.uint8
        return np.minimum(33 + qvArray, 126).tostring()

    def writeRecord(self, sequenceHeader, sequenceArray, qualityArray):
        """
        Sequence array is understood to be a string, or a numpy array
        containing strings which must be concatenated.  Quality array
        is a numpy array (dtype np.uint8) containing one number for
        each entry in sequenceArray.
        """
        assert len(sequenceArray) == len(qualityArray)
        self.file.write("@" + sequenceHeader + "\n")
        for chunk in chunks(OUTPUT_CHUNK_SIZE, sequenceArray):
            self.file.write("".join(chunk))
        self.file.write("\n")
        self.file.write("+")
        self.file.write("\n")
        for (seqChunk, qualityChunk) in itertools.izip(chunks(OUTPUT_CHUNK_SIZE, sequenceArray),
                                                       chunks(OUTPUT_CHUNK_SIZE, qualityArray)):
            lens = [len(s) for s in seqChunk]
            consensusQvs = np.repeat(qualityChunk, lens)
            self.file.write(self._qvArrayToString(consensusQvs))
        self.file.write("\n")

    def close(self):
        if not self.file.closed:
            self.file.close()

    def __del__(self):
        self.close()
