import string
from pypeline.common.makefile import \
    read_makefile, \
    IsUnsignedInt, \
    IsFloat, \
    IsStr, \
    IsBoolean, \
    And, \
    Or, \
    ValueIn, \
    IsNone, \
    StringIn, \
    ValueGE, \
    ValuesSubsetOf, \
    IsListOf

EXCLUDEBED = Or(IsListOf(IsStr), IsStr, IsNone, default=None)

def _alphanum_check(whitelist):
    description = "characters a-z, A-Z, 0-9%s allowed"
    description %= (", and %r" % whitelist,) if whitelist else ""

    whitelist += string.ascii_letters + string.digits

    return And(IsStr(),
               ValuesSubsetOf(whitelist, description=description))



_VALID_BED_NAME = _VALID_TARGET_NAME = \
    And(_alphanum_check(whitelist="._-"),
        ValueGE(2, key=len, description="at least two characters long"))

_VALIDATION_OPTIONS = {
    "BamPath": IsStr,
    "--MinMappingQuality": IsUnsignedInt(default=30),
    "--MinAlignmentLength": IsUnsignedInt(default=25),
    "--NoReadsChecked": IsUnsignedInt(default=10000)
}

_VALIDATION_GCCORRECT = {
    "Enabled": IsBoolean(default=False),
    "--NoRegions": Or(IsUnsignedInt, IsStr, default=200), ## Add a key to check all 
    "--ChromUsed": Or(IsStr,IsUnsignedInt, default="all"),   ## the is new
    "--MappaUniqueness": IsFloat(default=0.9)
}



_VALIDATION_NUCLEO = {
    "Enabled": IsBoolean(default=False),
    "Apply_GC_Correction": IsBoolean(default=True),
    "ExcludeBed": EXCLUDEBED,
    "--NucleosomeFlanks": IsUnsignedInt(default=25),
    "--NucleosomeSize": IsUnsignedInt(default=147),
    "--NucleosomeOffset": IsUnsignedInt(default=12)
}
_VALIDATION_METHYL = {
    "Enabled": IsBoolean(default=False),
    "ExcludeBed": EXCLUDEBED,
    "--ReadBases": IsUnsignedInt(default=15),
    "--SkipThreePrime": IsUnsignedInt(default=0),
    "--SkipFivePrime": IsUnsignedInt(default=0),
    "--Primes": StringIn(('five', 'three','both'), default='five')
}
_VALIDATION_PHASO = {
    "Enabled": IsBoolean(default=False),
    "ExcludeBed": EXCLUDEBED,
    "Apply_GC_Correction": IsBoolean(default=False),
    "--SubsetPileup": IsUnsignedInt(default=3),
    "--MaxRange": IsUnsignedInt(default=1000)
}
_VALIDATION_WRITEDEPTH = {
    "Enabled": IsBoolean(default=False),
    "ExcludeBed": EXCLUDEBED,
    "Apply_GC_Correction": IsBoolean(default=True)
}

_VALIDATION_TOOLS = {
    "BamInfo": _VALIDATION_OPTIONS,
    "GCcorrect": _VALIDATION_GCCORRECT,
    "NucleoMap": _VALIDATION_NUCLEO,
    "MethylMap": _VALIDATION_METHYL,
    "Phasogram": _VALIDATION_PHASO,
    "WriteDepth": _VALIDATION_WRITEDEPTH
}
# _VALIDATION_BED = {
#     "Path": IsStr
#     # "MakeMergePlot": IsBoolean(default=False)
# }

_VALIDATION = {
    "BamInputs": {  # BAMFILES:
        _VALID_TARGET_NAME: _VALIDATION_TOOLS
    },
    "Prefixes": {
        "--FastaPath": IsStr,
        "--MappabilityPath": Or(IsStr, IsNone)
    },
    "BedFiles": {
        'MappabilityFilter': IsBoolean(default=False),
        'MappabilityScore': IsFloat(default=0.9),
        _VALID_BED_NAME: IsStr
        ## _VALID_BED_NAME: _VALIDATION_BED
    }
}


def read_epiomix_makefile(argv):
    for f in argv:
        if isinstance(f, str) and f.endswith('.yaml'):
            yield read_makefile(f, _VALIDATION)
