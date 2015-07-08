#!/usr/bin/env python
import sys
import os
import time
import traceback
from optparse import OptionParser
from pbcore.io import GffReader, WriterBase

#
# (Ported from pbpy)
#

class VcfRecord:
    """Models a record in a VCF3.3 file."""
    def __init__(self):
        self.chrom = ''
        self.pos = 1
        self.id = '.'
        self.ref = ''
        self.alt = ''
        self.qual = -1.00
        self.filter = '0'
        self.info = {}

    @staticmethod
    def fromVariantGffRecord(gff):
        vcf = VcfRecord()
        vcf.chrom = gff.seqid
        vcf.id = '.'

        ref = gff.reference
        if ref is None:
            vcf.ref = "N"
        else:
            vcf.ref = ref
        vcf.qual = float(gff.confidence)
        vcf.put('NS', 1)
        vcf.put('DP', gff.coverage)

        feature = gff.type
        vcf.pos = gff.start
        if feature == 'insertion':
            vcf.alt = 'I%s' % gff.variantSeq.upper()
        elif feature == 'deletion':
            vcf.alt = 'D%s' % len(gff.reference)
        elif feature == 'substitution':
            vcf.alt = gff.variantSeq.upper()
        else:
            print >> sys.stderr, 'Unsupported feature %s found in GFF3 file.' % feature

        return vcf

    def put(self, key, value):
        self.info[key] = value

    @staticmethod
    def getHeader():
        return 'CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'

    def _getInfoString(self):
        return ';'.join(['%s=%s' % (k,v) \
            for k,v in self.info.iteritems()])

    def __str__(self):
        return '%s\t%d\t%s\t%s\t%s\t%.2f\t%s\t%s' % \
            (self.chrom, self.pos, self.id, self.ref, \
              self.alt, self.qual, self.filter, self._getInfoString())

class VcfWriter(WriterBase):
    """Outputs VCF (1000 Genomes Variant Call Format) 3.3 files"""
    def __init__(self, outfile):
        self._outfile = outfile
        self._start()

    def close(self):
        self._outfile.close()

    def flush(self):
        self._outfile.flush()

    def _start(self):
        self.writeMetaData('fileformat', 'VCFv3.3')

    def writeHeader(self):
        print >> self._outfile, '#%s' % VcfRecord.getHeader()

    def writeMetaData(self, key, value):
        print >> self._outfile, '##%s=%s' % (key, value)

    def writeRecord( self, record ):
        print >> self._outfile, str(record)

class GffToVcf(object):
    """Utility for converting variant GFF3 files to 1000 Genomes VCF"""
    def __init__(self, argv):
        self.__parseOptions(argv)

    def __parseOptions(self, argv):
        usage = 'Usage: %prog [--help] [options] variants.gff > variants.vcf'
        parser = OptionParser(usage=usage, description=GffToVcf.__doc__)
        parser.add_option('--globalReference', type="string",
                          help="Name of global reference to put in Meta field")
        parser.set_defaults(globalReference=None)

        self.options, self.args=parser.parse_args(argv)
        if len(self.args)!=2:
            parser.error('Expected 1 argument')

        self.gffFile = self.args[1]

    def _writeMetaData(self, writer):
        currentTime = time.localtime()
        cmdLine = os.path.basename(sys.argv[0]) + ' ' + ' '.join(sys.argv[1:])

        writer.writeMetaData('fileDate', '%d%d%d' % \
                             (currentTime[0], currentTime[1], currentTime[2]))
        writer.writeMetaData('source', cmdLine)
        if self.options.globalReference is not None:
            writer.writeMetaData('reference', self.options.globalReference)
        writer.writeMetaData('INFO', 'NS,1,Integer,"Number of Samples with Data"')
        writer.writeMetaData('INFO', 'DP,1,Integer,"Total Depth of Coverage"')

        writer.writeHeader()

    def run(self):
        with GffReader(self.gffFile) as reader, \
             VcfWriter(sys.stdout) as writer:
            self._writeMetaData(writer)
            for gff in reader:
                vcf = VcfRecord.fromVariantGffRecord(gff)
                writer.writeRecord(vcf)


if __name__ == '__main__':
    app = GffToVcf(sys.argv)
    sys.exit(app.run())
