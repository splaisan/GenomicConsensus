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

"""
I/O support for VCF output

The specification for the VCF format is available at
    http://www.1000genomes.org/wiki/Analysis/Variant%20Call%20Format/vcf-variant-call-format-version-41
"""

__all__ = [ "VcfWriter" ]

from .. import reference
from pbcore.io.base import WriterBase
import time

class VcfWriter(WriterBase):
    """
    A VCF file writer class specialized to GenomicConsensus's needs.
    Might be worth factoring it.
    """
    def __init__(self, f, sampleNames):
        super(VcfWriter, self).__init__(f)
        self.referenceFilename = options.referenceFilename
        self.sampleNames       = sampleNames
        self.writeDefaultHeaders()

    def writeDefaultHeaders(self):
        self.writeHeader("##fileformat=VCFv4.1")
        self.writeHeader("##source=GenomicConsensus")
        date = time.strftime("%Y%m%d", time.localtime())
        self.writeHeader("##fileDate=%s" % date)
        self.writeHeader("##reference=file://%s" % self.referenceFilename)

        for entry in self.referenceEntries:
            self.writeHeader("##contig=<ID=%s,length=%d,md5=%s>" %
                             (entry.name, entry.length, entry.md5sum))
        self.writeHeader("""##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">""")
        self.writeHeader("""##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">""")
        self.writeHeader("""##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">""")
        self.writeColumnNames()

    def writeHeader(self, headerLine):
        if not headerLine.startswith("##"):
            raise ValueError("VCF headers must start with ##")
        self.file.write("{0}\n".format(headerLine.rstrip()))

    def writeColumnNames(self):
        fixedFields = \
            ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
        genotypeFields = ["FORMAT"] + self.sampleNames
        self.file.write("#")
        self.file.write("\t".join(fixedFields + genotypeFields))
        self.file.write("\n")

    @staticmethod
    def _vcfEncodeVariants(refContig, refStart, varsAtPosition):
        # Identify (refPos, refAllele, altAlleles, genotypes), for
        # variant across samples, in the VCF convention, which
        # requires at least one base of context.
        #   - refPos is the 1-based start position of the padded reference allele
        #   - refAllele is the padded reference allele
        #   - altAlleles is the padded alternate allele sequences
        #   - genotypes is a list of the the encoded genotypes+info (GT:GQ:DP) for each sample
        refWindows = [vcfReferenceContext(v) for v in varsAtPosition]
        start = min(map(fst, refWindows))
        end   = max(map(snd, refWindows))
        refAllele = refContig[start:end]
        altAlleles = set(vcfAlternateAllele(refContig, start, end, var)
                         for var in varsAtPosition)

    def writeVariantsForPosition(self, refName, refStart, variantsAtPosition):
        # variantsAtPosition is a dict of
        #    Sample name -> Variant
        # if sample name is not present in the dict, it means
        # no variant was detected in the sample.

        # Assumptions:
        #   - All variants have same refStart
        #   - ?
        for sampleName in self.sampleNames:
            variantForSample = variantsAtPosition.get(sampleName)
            id_ = "."
            refAllele, altAlleles =

            self.file.write("%s\t%d\t"
                            % (refName, refStart+1))

def vcfReferenceContext(var):
    # Variants reported in VCF require a padding base... incredibly
    # obnoxious, this.  Normally the padding base is the base before
    # the event, but if the event occurs at the first base of the
    # genome, the padding base is the base after the event.
    # Context is returned as window [start, end) in 0-based half-open coordinates.
    if isinstance(Substitution):
        return var.refStart, var.refEnd
    elif var.refStart = 0:
        return var.refStart, var.refEnd + 1
    else:
        return var.refStart - 1, var.refEnd


def vcfAlternateAllele(refContig, start, end, var):
    # The padding potentially added to the reference context requires
    # recalculating the context of the variant sequence.
    padStart = var.refStart - start
    padEnd   = end - var.refEnd
    return refContig[start:var.refStart] + var.readSequence + refContig[var.refEnd:end]
