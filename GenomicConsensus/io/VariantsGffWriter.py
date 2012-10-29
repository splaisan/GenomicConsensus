
import time
from pbcore.io import GffWriter


class VariantsGffWriter(object):

    ONTOLOGY_URL = \
        "http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.12"

    def __init__(self, f, shellCommand, referenceEntries):
        self._gffWriter = GffWriter(f)
        self._gffWriter.writeMetaData("pacbio-variant-version", "1.4")
        self._gffWriter.writeMetaData("date", time.ctime())
        self._gffWriter.writeMetaData("feature-ontology", self.ONTOLOGY_URL)
        self._gffWriter.writeMetaData("source", "GenomicConsensus 0.3.0")
        self._gffWriter.writeMetaData("source-commandline",  shellCommand)

        # Reference groups.
        for entry in referenceEntries:
            self._gffWriter.writeMetaData("sequence-region", "%s 1 %d" \
                                          % (entry.header, entry.length))

    def writeRecord(self, gffRecord):
        self._gffWriter.writeRecord(gffRecord)
