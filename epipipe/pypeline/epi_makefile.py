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
    ValueGE, \
    ValuesSubsetOf, \
    IsListOf


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
    "--BamPrefix": Or(IsStr, ValueIn(('',)), default=''),
    "--LibraryConstruction": IsStr
}
_VALIDATION_GCCORRECT = {
    "Enabled": IsBoolean(default=True),
    "MapMinMaxReadLength": IsListOf(IsUnsignedInt),
    "--MappaUniqueness": IsFloat(default=0.9)
}
_VALIDATION_NUCLEO = {
    "Enabled": IsBoolean(default=True),
    "Apply_GC_Correction": IsBoolean(default=True),
    "--MinDepth": IsUnsignedInt(default=20),
    "--DequeLen": IsUnsignedInt(default=2000),
    "--NucleosomeFlanks": IsUnsignedInt(default=25),
    "--NucleosomeSize": IsUnsignedInt(default=147),
    "--NucleosomeOffset": IsUnsignedInt(default=12)
}
_VALIDATION_METHYL = {
    "Enabled": IsBoolean(default=True),
    "--ReadBases": IsUnsignedInt(default=6),
    "--SkipThreePrime": IsUnsignedInt(default=0),
    "--SkipFivePrime": IsUnsignedInt(default=0),
    "--Primes": IsStr
}
_VALIDATION_PHASO = {
    "Enabled": IsBoolean(default=True),
    "Apply_GC_Correction": IsBoolean(default=False),
    "--SubsetPileup": IsUnsignedInt(default=3),
    "--MaxRange": IsUnsignedInt(default=3000)
}
_VALIDATION_WRITEDEPTH = {
    "Enabled": IsBoolean(default=False),
    "Apply_GC_Correction": IsBoolean(default=True),
    "--DequeLength": IsUnsignedInt(default=500)
}

_VALIDATION_TOOLS = {
    "BamInfo": _VALIDATION_OPTIONS,
    "GCcorrect": _VALIDATION_GCCORRECT,
    "NucleoMap": _VALIDATION_NUCLEO,
    "MethylMap": _VALIDATION_METHYL,
    "Phasogram": _VALIDATION_PHASO,
    "WriteDepth": _VALIDATION_WRITEDEPTH
}

_VALIDATION_BED = {
    "Path": IsStr,
    "MakeMergePlot": IsBoolean(default=True),
    "NOSubFiles": IsUnsignedInt(default=3)
}

_VALIDATION = {
    "BamInputs": {  # BAMFILES:
        _VALID_TARGET_NAME: _VALIDATION_TOOLS
    },
    "Prefixes": {
        "--FastaPath": IsStr,
        "--FastaPrefix": Or(IsStr, ValueIn(('',)), default=''),
        "--MappabilityPath": Or(IsStr, ValueIn(('',)), default=''),
#        "--MappabilityPath": IsStr
    },
    "BedFiles": {
        'EnabledFilter': IsBoolean(default=True),
        'UniquenessFilter': IsFloat(default=0.9),
        _VALID_BED_NAME: _VALIDATION_BED
    }
}


def read_epiomix_makefile(argv):
    for f in argv:
        if isinstance(f, str) and f.endswith('.yaml'):
            yield read_makefile(f, _VALIDATION)
