IsStr, IsDictOf, Or = object(), object(), object()
BAMNAME = object()
_VALID_TARGET_NAME = object()
_VALDIATION_OPTIONS = object()
_VALID_PREFIX_NAME, _VALID_PREFIX_PATH = object(), object()
_VALID_MAPPA_NAME, _VALID_BED_NAME = object(), object()
IsUnsignedInt, IsBoolean, IsUnsignedInt = object(), object(), object()
ValuesSubsetOf, IsNone = object(), object()
And, ValueGE = object(), object()
_alphanum_check = object()

_VALIDATION_OPTIONS = {
    "--MinMappingQuality": IsUnsignedInt(default=30),
    "--BamChromType": IsStr(default=''),
    "--MapReadLengthMin": IsUnsignedInt,
    "--MapReadLengthMax": IsUnsignedInt,
    "--CorrectGCbias": IsBoolean(default=True),
    "--LibraryConstruction": IsStr
}
_VALIDATION_NUCLEO = {
    "Enabled": IsBoolean(default=True),
    "--MinDepth": IsUnsignedInt(default=20),
    "--DequeLen": IsUnsignedInt(default=2000)
}
_VALIDATION_METHYL = {
    "Enabled": IsBoolean(default=True),
    "--ReadBases": IsUnsignedInt(default=6),
    "--SkipBases": IsUnsignedInt(default=0)
}
_VALIDATION_PHASO = {
    "Enabled": IsBoolean(default=True),
    "--SubsetPileup": IsUnsignedInt(default=3),
    "--MaxRange": IsUnsignedInt(default=3000)
}
_VALIDATION_WRITEDEPTH = {
    "Enabled": IsBoolean(default=False)
}

_VALID_TARGET_NAME = \
    And(_alphanum_check(whitelist="._-"),
        ValueGE(2, key=len, description="at least two characters long"))


_VALIDATION_TOOLS = {
    "BAMNAME": IsStr,
    "Options": _VALIDATION_OPTIONS,
    "NucleoMap": _VALIDATION_NUCLEO,
    "MethylMap": _VALIDATION_METHYL,
    "Phasogram": _VALIDATION_PHASO,
    "WriteDepth": _VALIDATION_WRITEDEPTH
}

_VALIDATION = {
    _VALID_TARGET_NAME: {  # BAMFILES:
        _VALID_TARGET_NAME: _VALIDATION_TOOLS
    },
    "Prefixes": {
        _VALID_PREFIX_NAME: {
            "FastaPath": _VALID_PREFIX_PATH,
            "ChromType": IsStr
        },
        _VALID_MAPPA_NAME: {
            "MappabilityPath": IsStr
        }
    },
    "BedFiles": {
        _VALID_BED_NAME: {
            "BedPath": Or(IsStr, IsDictOf(IsStr, IsStr))
        }
    }
}
