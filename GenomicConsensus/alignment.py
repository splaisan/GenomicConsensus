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
