#!/usr/bin/env python
from __future__ import print_function
import sys
import datetime

_TEMPLATE_EPIPALEOMIX = \
"""# -*- mode: Yaml; -*-
# Timestamp: %s
# Please respect indentation (with spaces), and pay attention to colons and hyphens.
# Hash-commented lines are ignored.
# epiPALEOMIX requires as a minimum:
#   An indexed Reference genome
#   An indexed BAM file
#   A bedfile with region(s) of interest


# Only Names To Modify:
#   Bednames (i.e NameofBed)
#   Bamnames (i.e. BAMName1)

# incase of several bamfiles:
    # Add another 'BAMName2:' with a meaningful BAM file name.
    # Duplicate the info of the analyses in BAMName1 and change the required arguments (e.g. paths, and options)

Prefixes:
    --FastaPath: path/to/referenceFastafile           # REQUIRED
    
    # path to mappability file, leave empty if none
    --MappabilityPath:   # pathto/mappabilityfile

# At least one bed file is required
BedFiles:

    # "MappabilityFilter" {True, False, default=False}
        # If True, only uniquelymappable regions are analyzed. Requires valied MappabilityPath
        # If False, all bedregions will be used
    MappabilityFilter: False
    
    # "MappabilityScore" {a float between 0-1}
        # The uniqueness filter used.
    MappabilityScore: 0.9  # filtering bed regions with low uniqueness
    
    # For Each Bedfile Provided Set A Meaningful Bedname
    Nameofbed: path/to/bedfile
    # Nameofbed2: path/to/bedfile


BamInputs:
    # For Each BamFile Provided Set A Meaningful Bamname
    BAMName1:
        # "BamInfo" Contain general options about the Bamfile
        BamInfo:
            # Path to the bamfile
            BamPath: path/to/bamfile  # REQUIRED

            # Minimum Mapping Quality Read
            --MinMappingQuality: 25

            # Minimum Aligned Read Length
            --MinAlignmentLength: 25

            # Number of reads checked to identify minimum and maximum readlength. 10000 is enough
            --NoReadsChecked: 10000

        # "GCcorrection" contain options regarding the GC-correction model
        GCcorrect:
            # "Enabled" {True, False, default=False}. Requires a valid mappability file
            Enabled: False

            # The reference uniqueness minimum. It is not used if no mappability region is passed
            --MappaUniqueness: 0.9
            
            # If a specific chromsome should be used to create the gc-correction model
            # Default uses regions from entire genome
            --ChromUsed: all

            # Number of Regions with sufficiently high uniqueness to be investigated. 
            # Default is 200, that is sufficent for a 10x genome.
            # In case too few regions are being analysed, epiPALEOMIX will raise an except.
            # "--NoRegions" takes an integer or "all". Note if "all" can be timeconsuming.
            --NoRegions: 200

        # "NucleoMap" contain options regarding the nucleosome calling method
        NucleoMap:
            # "Enabled" {True, False, default=False}
            Enabled: False
           
            # "Apply_GC_Correction" {True, False, default=False}.
            # Whether the GC-correction model should be used when calling Nucleosomes
            Apply_GC_Correction: False  # GCcorrect Must be Enabled: True if applied

            # In case some of the Bedfiles should not be analyzed be NucleoMap
            # Either ExcludeBed: bed1 OR ExcludeBed: [bed1, bed2]
            ExcludeBed:  
            
            # "--NucleosomeSize" reflects number of nucleotides a nucleosome cover
            --NucleosomeSize: 147
            # "--NucleosomeFlanks" reflects number of nucleotides used to calculated the neighbouring coverage
            --NucleosomeFlanks: 25
            # "--NucleosomeOffset" reflects number of nucleotides between Flanks and Size
            --NucleosomeOffset: 12

        # "MethylMap" contain options regarding the methylation mapping tool
        MethylMap:
            # "Enabled" {True, False, default=False}
            Enabled: False

            # In case some of the Bedfiles should not be analyzed be NucleoMap
            # Either ExcludeBed: bed1 OR ExcludeBed: [bed1, bed2]
            ExcludeBed:  

            # "--Primes" regards what prime of the reads should be investigated.
            # For doublestranded library protocols  'five' (i.e. checking five prime end of reads only) should be used. 
            # For single stranded library protocols 'both' (i.e. checking both primes of the reads)
            # 'three' (i.e. checking only three prime end of reads) is also available but rarely used.
            --Primes: five

            # Number of bases being analyzed.
            --ReadBases: 15

            # In cases nucleotides  from the three or five prime should be excluded.
            # The number of nucleotides removed from the primes are subtracted from --ReadBases.
            --SkipThreePrime: 0
            --SkipFivePrime: 0

        # "Phasogram" contain options regarding the Phasogram tool
        Phasogram:
            # "Enabled" {True, False, default=False}
            Enabled: False

            # In case some of the Bedfiles should not be analyzed be NucleoMap
            # Either ExcludeBed: bed1 OR ExcludeBed: [bed1, bed2]
            ExcludeBed:  

            # "Apply_GC_Correction" {True, False, default=False}.
            # Whether the GC-correction model should be used when doing the phasogram analysis
            Apply_GC_Correction: False  # GCcorrect Must be "Enabled: True" if this is applied

            # "--SubsetPileup" reflects the number of reads required to share 5' startposition. 
            --SubsetPileup: 3

            # "--MaxRange" reflects the distance in which the phasogram should be calculated
            --MaxRange: 1000

        # "WriteDepth" contain options regarding the Pileup depth and running-nucleosome score tool
        WriteDepth:
            # "Enabled" {True, False, default=False}
            Enabled: False

            # In case some of the Bedfiles should not be analyzed be NucleoMap
            # Either ExcludeBed: bed1 OR ExcludeBed: [bed1, bed2]
            ExcludeBed:  

            # "Apply_GC_Correction" {True, False, default=False}.
            # Whether the GC-correction model should be used when doing the phasogram analysis
            Apply_GC_Correction: False  # GCcorrect Must be "Enabled: True" if this is applied
"""
_TEMPLATE_EPIPALEOMIX_SIMPLE = \
"""# -*- mode: Yaml; -*-
# Timestamp: %s
# Please respect indentation (with spaces), and pay attention to colons and hyphens.
# Hash-commented lines are ignored.
Prefixes:
    --FastaPath: path/to/referenceFastafile           # REQUIRED
    --MappabilityPath: # pathto/mappabilityfile
BedFiles:
    MappabilityFilter: False
    MappabilityScore: 0.9  # filtering bed regions with low uniqueness
    Nameofbed: path/to/bedfile
    # Nameofbed2: path/to/bedfile
BamInputs:
    BAMName1:
        BamInfo:
            BamPath: path/to/bamfile  # REQUIRED
        GCcorrect:
            Enabled: False
        NucleoMap:
            Enabled: False
            Apply_GC_Correction: False  # GCcorrect Must be Enabled: True if applied
        MethylMap:
            Enabled: False
        Phasogram:
            Enabled: False
            Apply_GC_Correction: False  # GCcorrect Must be Enabled: True if applied
        WriteDepth:
            Enabled: False
            Apply_GC_Correction: False  # GCcorrect Must be Enabled: True if applied
"""


def main(argv):
    timestamp = datetime.datetime.now().isoformat()
    try:
        if argv[0].lower() == 'simple':
            template = _TEMPLATE_EPIPALEOMIX_SIMPLE % (timestamp,)
            print(template, file=sys.stdout)
            return 0
    except IndexError:
        pass
    template = _TEMPLATE_EPIPALEOMIX % (timestamp,)
    print(template, file=sys.stdout)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
