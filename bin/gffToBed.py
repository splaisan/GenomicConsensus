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

class BedRecord:
    """Models a record in a BED file format"""
    def __init__(self):
        self.chrom=''
        self.chromStart = 0
        self.chromEnd = 0
        self.name = ''
        self.score = -1.00
        self.strand = '+'

    def __str__(self):
        return '%s\t%d\t%d\t%s\t%.3f\t%s' % \
            (self.chrom, self.chromStart, self.chromEnd, self.name, \
              self.score, self.strand)

class CoverageBedRecord(BedRecord):
    @staticmethod
    def fromAlignmentSummaryGffRecord(gff):
        bed = CoverageBedRecord()
        bed.chrom = gff.seqid
        bed.chromStart = gff.start - 1
        bed.chromEnd = gff.end
        bed.name = 'meanCov'
        bed.score = float(gff.cov2.split(',')[0])
        bed.strand = gff.strand
        return bed

class VariantsBedRecord(BedRecord):
    @staticmethod
    def fromVariantGffRecord(gff):
        bed = VariantsBedRecord()
        bed.chrom = gff.seqid
        bed.chromStart = gff.start - 1
        bed.score = float(gff.confidence)
        bed.strand = gff.strand

        feature = gff.type
        #GFF3 coordinates are 1-based and inclusive
        #BED coordinates are 0-based and exclusive
        if feature == 'insertion':
            bed.chromEnd = bed.chromStart + 1
            bed.name = '%d_%dins%s' % (bed.chromStart + 1,
                                       bed.chromEnd + 1,
                                       gff.variantSeq)
        elif feature == 'deletion':
            featureLen = len(gff.reference)
            bed.chromEnd = bed.chromStart + featureLen
            if featureLen == 1:
                bed.name = "%ddel" % (bed.chromStart + 1)
            else:
                bed.name = '%d_%ddel' % (bed.chromStart + 1, bed.chromEnd)
        elif feature == 'substitution':
            bed.chromEnd = bed.chromStart + 1
            bed.name = '%d%s>%s' % (bed.chromStart + 1,
                                    gff.reference,
                                    gff.variantSeq)
        else:
            print >> sys.stderr, 'Unsupported feature %s found in GFF3 file.' % feature

        return bed

class BedWriter(WriterBase):
    """Outputs BED annotation track file"""
    def __init__(self, outfile):
        self._outfile = outfile

    def close(self):
        self._outfile.close()

    def flush(self):
        self._outfile.flush()

    def writeHeader(self, name, description, useScore):
        print >> self._outfile, 'track name=%s description="%s" useScore=%d' \
            % (name, description, useScore)

    def writeRecord(self, record):
        print >> self._outfile, str(record)

class GffToBed:
    """
    Utility for converting GFF3 to BED format. Currently supports
    regional coverage or variant .bed output.
    """
    def __init__(self, argv):
        self.__parseOptions(argv)

    def __parseOptions(self, argv):
        usage = 'Usage: %prog [--help] [options] purpose[coverage|variants] input.gff > output.bed'
        parser = OptionParser(usage=usage, description=GffToBed.__doc__)

        parser.add_option('--name', type="string",
                          help="track name to display in header")
        parser.add_option('--description', type="string",
                          help="track description to display in header")
        parser.add_option('--useScore', type="int", default=0,
                          help="whether or not to use score for feature display")

        self.options, self.args=parser.parse_args(argv)
        if len(self.args)!=3:
            parser.error('Expected 2 argument')

        self.purpose = self.args[1]

        if self.purpose not in [ "variants", "coverage" ]:
            print >> sys.stderr, \
                "Purpose %s not supported. Must be one of: [variants|coverage]" % (self.purpose)
            sys.exit(-1)

        self.gffFile = self.args[2]

    def run(self):
        with GffReader(self.gffFile) as reader, \
             BedWriter(sys.stdout)   as writer:

            writer.writeHeader(self.options.name,
                               self.options.description,
                               self.options.useScore)
            for gff in reader:
                if self.purpose == 'coverage':
                    bedRecord = CoverageBedRecord.fromAlignmentSummaryGffRecord(gff)
                else:
                    bedRecord = VariantsBedRecord.fromVariantGffRecord(gff)
                writer.writeRecord(bedRecord)


if __name__ == '__main__':
    app = GffToBed(sys.argv)
    sys.exit(app.run())
