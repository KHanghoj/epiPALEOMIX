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


# Only Names To Modify. (Do not use underscores):
#   Bednames (i.e NameofBed)
#   Bamnames (i.e. BAMName1)

# in-case of several bamfiles:
    # Add another 'BAMName2:' with a meaningful BAM file name.
    # Keep the second bamfile on the same indentation level as BAMName1

Prefixes:
    --FastaPath: path/to/referenceFastafile       # REQUIRED
    
    # path to mappability file, leave empty if none.
    # GCcorrection cannot be conducted with out a mappability file.
    --MappabilityPath:   # pathto/mappabilityfile

# At least one bed file is required
BedFiles:

    # "MappabilityFilter" {True, False, default=False}
        # If True, only uniquely-mappable regions are analyzed. Requires valid MappabilityPath
        # If False, all bedregions will be used
    MappabilityFilter: False
    
    # "MappabilityScore" {a float between 0-1}
        # The uniqueness score used to filter.
    MappabilityScore: 0.9  # filtering bed regions with low uniqueness
    
    # For Each Bedfile Provided Set A Meaningful Bedname
    Nameofbed: path/to/bedfile                    # REQUIRES at least one bedfile
    # Nameofbed2: path/to/bedfile


BamInputs:
    # For Each BamFile Provided Set A Meaningful Bamname
    BAMName1:
        # "BamInfo" Contain general options about the Bamfile
        BamInfo:
            # Path to the bamfile
            BamPath: path/to/bamfile              # REQUIRED

            # Minimum Mapping Quality Read
            --MinMappingQuality: 25

            # Minimum Aligned Read Length
            --MinAlignmentLength: 25

            # Number of reads checked to identify minimum and maximum readlength, roughly. 5000 is enough
            --NoReadsChecked: 5000

        # "GCcorrection" contain options regarding the GC-correction model
        GCcorrect:
            # "Enabled" {True, False, default=False}. 
            #Requires a valid mappability file '--mappabilityPath' in Prefix
            Enabled: False

            # The reference uniqueness minimum. It is not used if no mappability region is passed
            --MappaUniqueness: 0.9
            
            # If a specific chromsome should be used to create the gc-correction model
            # Default uses regions from entire genome
            --ChromUsed: all

            # Number of Regions with sufficiently high uniqueness to be investigated. 
            # Default is 200, that is sufficent for a 5x genome.
            # In case too few regions are being analyzed, epiPALEOMIX will raise an exception.
            # "--NoRegions" takes an integer or "all". Note that using  "all" will be timeconsuming.
            --NoRegions: 200

        # "NucleoMap" contains options for nucleosome calling
        NucleoMap:
            # "Enabled" {True, False, default=False}
            Enabled: False
           
            # "Apply_GC_Correction" {True, False, default=False}.
            # Whether the GC-correction model should be used when calling Nucleosomes
            Apply_GC_Correction: False  # GCcorrect Must be Enabled: True if applied

            # In case some of the Bedfiles should not be analyzed by NucleoMap
            # Either ExcludeBed: bed1 OR ExcludeBed: [bed1, bed2]
            ExcludeBed:  
            
            # "--NucleosomeSize" reflects the number of nucleotides a nucleosome cover
            --NucleosomeSize: 147
            # "--NucleosomeFlanks" reflects the number of nucleotides used to calculated the neighbouring coverage
            --NucleosomeFlanks: 25
            # "--NucleosomeOffset" reflects the number of nucleotides between Flanks and Size
            --NucleosomeOffset: 12

        # "MethylMap" contains options the methylation mapping tool
        MethylMap:
            # "Enabled" {True, False, default=False}
            Enabled: False

            # In case some of the Bedfiles should not be analyzed be MethylMap
            # Either ExcludeBed: bed1 OR ExcludeBed: [bed1, bed2]
            ExcludeBed:  

            # "--Primes" indicates which read termini (5-prime and/or 3-prime) should be investigated:
            # For doublestranded library protocols  'five' (i.e. 5-prime) should be used. 
            # For single stranded library protocols 'both' i.e. both 5-prime and 3-prime) should be used
            # 'three' (i.e. checking only three prime end of reads) is also available but rarely used.
            --Primes: five

            # Number of bases being analyzed.
            --ReadBases: 15

            # In cases nucleotides  from the three or five prime should be excluded.
            # The number of nucleotides removed from the primes are subtracted from --ReadBases, and can, thus, not be greater than ReadBases.
            --SkipThreePrime: 0
            --SkipFivePrime: 0

        # "Phasogram" contains options for the Phasogram tool
        Phasogram:
            # "Enabled" {True, False, default=False}
            Enabled: False

            # In case some of the Bedfiles should not be analyzed by Phasogram
            # Either ExcludeBed: bed1 OR ExcludeBed: [bed1, bed2]
            ExcludeBed:  

            # "Apply_GC_Correction" {True, False, default=False}.
            # Whether the GC-correction model should be used when doing the phasogram analysis
            Apply_GC_Correction: False  # GCcorrect Must be "Enabled: True" if this is applied

            # "--SubsetPileup" reflects the minimal number of reads required to share 5' start positions.
            --SubsetPileup: 1

            # "--MaxRange" reflects the distance (in bp) within which the phasogram should be calculated
            --MaxRange: 1000

        # "WriteDepth" contains options for the Pileup depth and running-nucleosome score tools
        WriteDepth:
            # "Enabled" {True, False, default=False}
            Enabled: False

            # In case some of the Bedfiles should not be analyzed by WriteDepth
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
    --FastaPath: path/to/referenceFastafile       # REQUIRED
    --MappabilityPath: # pathto/mappabilityfile
BedFiles:
    MappabilityFilter: False
    MappabilityScore: 0.9                         # filtering bed regions with low uniqueness
    Nameofbed: path/to/bedfile                    # REQUIRED
    # Nameofbed2: path/to/bedfile
BamInputs:
    BAMName1:
        BamInfo:
            BamPath: path/to/bamfile              # REQUIRED
        GCcorrect:
            Enabled: False
        NucleoMap:
            Enabled: False
            Apply_GC_Correction: False            # GCcorrect Must be Enabled if applied
        MethylMap:
            Enabled: False
        Phasogram:
            Enabled: False
            Apply_GC_Correction: False            # GCcorrect Must be Enabled if applied
        WriteDepth:
            Enabled: False
            Apply_GC_Correction: False            # GCcorrect Must be Enabled if applied
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
