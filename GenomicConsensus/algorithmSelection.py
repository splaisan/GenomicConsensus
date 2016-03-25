#!/usr/bin/env python
#################################################################################
# Copyright (c) 2011-2016, Pacific Biosciences of California, Inc.
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

def bestAlgorithm(sequencingChemistries, logger=None):
    """
    Identify the (de novo) consensus algorithm we expect to deliver
    the best results, given the sequencing chemistries represented in
    an alignment file.

    We key off the sequencing chemistries as follows:

    - Just RS chemistry data?  Then use quiver (at least for now, until
      we get arrow > quiver on P6-C4)
    - S/P1-C1/beta chemistry present?  Then we have to fall back to POA.
    - Else (either all Sequel data post beta, or a mix of Sequel and
      RS data), use arrow.

    Note that the handling/rejection of chemistry mixtures (including
    mixtures of Sequel and RS data) is left to the algorithm itself.
    """
    if len(sequencingChemistries) == 0:
        raise ValueError("sequencingChemistries must be nonempty list or set")
    chems = set(sequencingChemistries)
    anyUnknown    = "unknown" in chems
    anySequelBeta = "S/P1-C1/beta" in chems
    allRS         = all(not(chem.startswith("S/")) for chem in chems) and (not anyUnknown)

    if anyUnknown:
        if logger is not None:
            logger.warning("Unidentifiable sequencing chemistry present in dataset; falling back to " +
                           "nonparametric but suboptimal 'POA' algorithm.  Check if your SMRTanalysis " +
                           "installation is out-of-date.")
        return "poa"
    elif anySequelBeta:
        return "poa"
    elif allRS:
        return "quiver"
    else:
        return "arrow"
