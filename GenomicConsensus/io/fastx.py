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
