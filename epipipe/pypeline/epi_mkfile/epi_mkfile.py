#!/usr/bin/python
from __future__ import print_function
import sys
import datetime

_TEMPLATE_EPIPALEOMIX = \
"""# -*- mode: Yaml; -*-
# Timestamp: %s

#DOT NOT CHANGES OPTIONS UNLESS OTHERWISE NOTED
#DOT NOT CHANGES THE NAMES ON THIS INDENTATION LEVEL
Prefixes:
    "--FastaPath": path/to/referenceFastafile  # REQUIRED
    "--FastaPrefix": ''  # if annotaed by digits only ('1')  leave it '' else add the prefix i.e. 'chr'
    "--MappabilityPath": path/to/Mappabilityfile # leave empty ('') if no mappability file
BedFiles:  # AT LEAST ONE BEDFILE IS REQUIRED
    EnabledFilter: False  # if True, bed coordinates not overlapping a mappability region with at at least 0.9 in uniqueness (user definable)
    UniquenessFilter: 0.9  # filtering bed regions with low uniqueness
    Nameofbed: # SET TO A MEANINGFUL BEDNAME
        Path: path/to/bedfile
        MakeMergePlot: False  # i.e. if bedfile covers chr1  cannot merge data
        NOSubFiles: 5 # splits the bedfile inputs into subfiles to take advantage of all cpu's
    # Nameofbed1: # SET TO A MEANINGFUL BEDNAME
    #     Path: path/to/bedfile
    #     MakeMergePlot: False  # i.e. if bedfile covers chr1  cannot merge data
    #     NOSubFiles: 5  # splits the bedfile inputs into subfiles to take advantage of all cpu's
BamInputs:
    BAMName1:  # SET TO A MEANINGFUL BAMNAME 
        # DOT NOT CHANGES THE VARIABLE NAMES ON THIS INDENTATION LEVEL
        BamInfo:
            "BamPath": path/to/bamfile  # REQUIRED
            "--MinMappingQuality": 30
            "--BamPrefix": ''  # if annotaed by number only '1'  leave it '' else add the prefix i.e. 'chr' 
            "--LibraryConstruction": DS  # DS if double strand, SS if single strand
        GCcorrect:
            "Enabled": True
            # mininimum and maximum length of used reads (roughly)
            "MapMinMaxReadLength": [30, 75] 
            # The reference uniqueness minimum. It is not used if no mappability region is passed
            "--MappaUniqueness": 0.9
        NucleoMap:   # The nucleosome calling tool
            "Enabled": True
            "Apply_GC_Correction": True  # GCcorrect Enabled must be True
            "--MinDepth": 20  # the minimum depth for a called nucleosome. suggested defualy is X coverage
            "--DequeLen": 2000
            "--NucleosomeSize": 147
            "--NucleosomeFlanks": 25
            "--NucleosomeOffset": 12
        MethylMap:   # The Methylation mapping tool
            "Enabled": True
            "--ReadBases": 6
            "--Primes": five  # choose from {both, five, three}. Normally 5' (five) for DS and both for SS
            "--SkipThreePrime": 0
            "--SkipFivePrime": 0
        Phasogram:   # THe phasogram analysis 
            "Enabled": True
            "Apply_GC_Correction": False
            "--SubsetPileup": 3
            "--MaxRange": 1000
        WriteDepth:
            "Enabled": False
            "Apply_GC_Correction": True
            "--DequeLength": 500
    # BAMName2:  # SET TO A MEANINGFUL BAMNAME 
    #     # DOT NOT CHANGES THE VARIABLE NAMES ON THIS INDENTATION LEVEL
    #     BamInfo:
    #         # DOT NOT CHANGES THE VARIABLE NAMES ON THIS INDENTATION LEVEL
    #         "BamPath": path/to/bamfile
    #         "--MinMappingQuality": 30
    #         "--BamPrefix": ''  # if annotaed by number only '1'  leave it '' else add the prefix i.e. 'chr' 
    #         "--LibraryConstruction": DS  # DS if double strand, SS if single strand
    #     GCcorrect:
    #         "Enabled": True
    #         # mininimum and maximum length of used reads (roughly)
    #         "MapMinMaxReadLength": [30, 75] 
    #         "--MappaUniqueness": 0.9  # The reference uniqueness minimum
    #     NucleoMap:   # The nucleosome calling tool
    #         "Enabled": True
    #         "Apply_GC_Correction": True
    #         "--MinDepth": 20  # the minimum depth for a called nucleosome. suggested defualy is X coverage
    #         "--DequeLen": 2000
    #         "--NucleosomeSize": 147
    #         "--NucleosomeFlanks": 25
    #         "--NucleosomeOffset": 12
    #     MethylMap:
    #         "Enabled": True
    #         "--ReadBases": 6
    #         "--Primes": five  # choose from {both, five, three}. Normally 5' for DS and both for SS
    #         "--SkipThreePrime": 0
    #         "--SkipFivePrime": 0
    #     Phasogram:
    #         "Enabled": True
    #         "Apply_GC_Correction": False
    #         "--SubsetPileup": 3
    #         "--MaxRange": 1000
    #     WriteDepth:
    #         "Enabled": False
    #         "Apply_GC_Correction": True
    #         "--DequeLength": 500
"""


def main(argv):
    timestamp = datetime.datetime.now().isoformat()
    template = _TEMPLATE_EPIPALEOMIX % (timestamp,)
    print(template, file=sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
