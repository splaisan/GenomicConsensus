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

# Container for alignment-related data structures.
import heapq
from collections import defaultdict

class AlignmentColumn(object):
    """
    An AlignmentColumn is a frequency tables of read snippets seen at
    a given reference position, allowing for quickly identifying the
    most common read snippet, and for identifying the confidence of
    our base call for this position.

    Unlike `collections.Counter`, this class has well-defined
    tie-breaking behavior (largest snippet, in lexicrographical
    order) for the "most frequent snippet."
    """

    def __init__(self, referenceId, referencePos, referenceBase):
        self.referenceId = referenceId
        self.referencePos = referencePos
        self.referenceBase = referenceBase
        self._snippetsAndFrequencies = defaultdict(int)
        # Ordered Snippets And Frequencies :-)
        self._osaf = None

    @staticmethod
    def _itemComparisonKey(item):
        """
        Comparisons are done first on frequency, and then
        on the snippet to break ties.  The lexicographically
        larger snippet wins in ties.
        """
        return (item[1], item[0])

    def __repr__(self):
        response = "AlignmentColumn [reference=%s]:\n" % self.referenceBase
        for snippet, freq in self.orderedSnippetsAndFrequencies:
            response += "%8s : %d" % (snippet, freq) + "\n"
        return response

    def addReadSnippet(self, snippet):
        self._snippetsAndFrequencies[snippet] += 1

    def coverage(self):
        return sum(self._snippetsAndFrequencies.values())

    def mostFrequentItem(self):
        return heapq.nlargest(1, self._snippetsAndFrequencies.iteritems(),
                              key=self._itemComparisonKey)[0]

    @property
    def orderedSnippetsAndFrequencies(self):
        # initialize lazily
        if not self._osaf:
            self._osaf = sorted(self._snippetsAndFrequencies.iteritems(),
                          key=self._itemComparisonKey)
            self._osaf.reverse()

        return self._osaf

    def confidence(self):
        pass

    def consensus(self):
        pass
