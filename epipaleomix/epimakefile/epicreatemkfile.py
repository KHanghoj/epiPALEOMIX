#!/usr/bin/python
from __future__ import print_function
import sys
import datetime

_TEMPLATE_EPIPALEOMIX = \
"""# -*- mode: Yaml; -*-
# Timestamp: %s
#DOT NOT CHANGES THE NAMES ON THIS INDENTATION LEVEL
Prefixes:
    --FastaPath: path/to/referenceFastafile  # REQUIRED
    --MappabilityPath:  # path to mappability file, leave empty if none
BedFiles:  # AT LEAST ONE BEDFILE IS REQUIRED
    EnabledFilter: False  # if True, bed coordinates not overlapping a mappability region with at at least 0.9 in uniqueness (user definable)
    UniquenessFilter: 0.9  # filtering bed regions with low uniqueness
    Nameofbed: path/to/bedfile 
    # Nameofbed1: # SET TO A MEANINGFUL BEDNAME
    #     Path: path/to/bedfile
BamInputs:
    BAMName1:  # SET TO A MEANINGFUL BAMNAME 
        # DOT NOT CHANGES THE VARIABLE NAMES ON THIS INDENTATION LEVEL
        BamInfo:
            BamPath: path/to/bamfile  # REQUIRED
            --MinMappingQuality: 30
        GCcorrect:
            Enabled: True
            # mininimum and maximum length of used reads (roughly)
            # ChromUsed: 1
            --NoRegions: 200
            # The reference uniqueness minimum. It is not used if no mappability region is passed
            --MappaUniqueness: 0.9
        NucleoMap:   # The nucleosome calling tool
            Enabled: True
            Apply_GC_Correction: True  # GCcorrect Enabled must be True
            --NucleosomeSize: 147
            --NucleosomeFlanks: 25
            --NucleosomeOffset: 12
        MethylMap:   # The Methylation mapping tool
            Enabled: True
            --ReadBases: 15
            --Primes: five  # choose from {both, five, three}. Normally 5' (five) for DS and 5' and 3' (both) for SS
            --SkipThreePrime: 0
            --SkipFivePrime: 0
        Phasogram:   # The phasogram analysis 
            Enabled: True
            Apply_GC_Correction: False
            --SubsetPileup: 3
            --MaxRange: 1000
        WriteDepth:
            Enabled: False
            Apply_GC_Correction: True
            --DequeLength: 500
    # BAMName2:  # SET TO A MEANINGFUL BAMNAME 
    #     # DOT NOT CHANGES THE VARIABLE NAMES ON THIS INDENTATION LEVEL
    #     BamInfo:
    #         BamPath: path/to/bamfile  # REQUIRED
    #         --MinMappingQuality: 30
    #     GCcorrect:
    #         Enabled: True
    #         # mininimum and maximum length of used reads (roughly)
    #         ChromUsed: [1,1]
    #         # The reference uniqueness minimum. It is not used if no mappability region is passed
    #         --MappaUniqueness: 0.9
    #     NucleoMap:   # The nucleosome calling tool
    #         Enabled: True
    #         Apply_GC_Correction: True  # GCcorrect Enabled must be True
    #         --NucleosomeSize: 147
    #         --NucleosomeFlanks: 25
    #         --NucleosomeOffset: 12
    #     MethylMap:   # The Methylation mapping tool
    #         Enabled: True
    #         --ReadBases: 15
    #         --Primes: five  # choose from {both, five, three}. Normally 5' (five) for DS and 5' and 3' (both) for SS
    #         --SkipThreePrime: 0
    #         --SkipFivePrime: 0
    #     Phasogram:   # The phasogram analysis 
    #         Enabled: True
    #         Apply_GC_Correction: False
    #         --SubsetPileup: 3
    #         --MaxRange: 1000
    #     WriteDepth:
    #         Enabled: False
    #         Apply_GC_Correction: True
    #         --DequeLength: 500
"""


def main(argv):
    timestamp = datetime.datetime.now().isoformat()
    template = _TEMPLATE_EPIPALEOMIX % (timestamp,)
    print(template, file=sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
