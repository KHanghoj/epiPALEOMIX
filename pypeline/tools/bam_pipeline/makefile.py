#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import os
import sys
import copy
import glob
import types
import string
import itertools
import collections

import pysam

import pypeline.tools.bam_pipeline.paths as paths
from pypeline.common.utilities import fill_dict
from pypeline.common.fileutils import missing_files
from pypeline.common.makefile import \
    MakefileError, \
    REQUIRED_VALUE, \
    WithoutDefaults, \
    read_makefile, \
    IsInt, \
    IsUnsignedInt, \
    IsFloat, \
    IsStr, \
    IsNone, \
    IsBoolean, \
    And, \
    Or, \
    Not, \
    ValueGE, \
    ValueIn, \
    ValuesIntersect, \
    ValuesSubsetOf, \
    StringIn, \
    StringStartsWith, \
    IsListOf, \
    IsDictOf
from pypeline.common.console import \
    print_info, \
    print_warn

import pypeline.nodes.bwa as bwa
import pypeline.common.versions as versions
import pypeline.common.sequences as sequences


_READ_TYPES = set(("Single", "Collapsed", "CollapsedTruncated", "Paired"))


def read_makefiles(config, filenames):
    makefiles = []
    for filename in filenames:
        makefile = read_makefile(filename, _VALIDATION)
        makefile = _mangle_makefile(makefile)

        makefiles.append(makefile)

    return _validate_makefiles(config, makefiles)


def _alphanum_check(whitelist):
    description = "characters a-z, A-Z, 0-9%s allowed"
    description %= (", and %r" % whitelist,) if whitelist else ""

    whitelist += string.ascii_letters + string.digits

    return And(IsStr(),
               ValuesSubsetOf(whitelist, description=description))


# Valid names for prefixes
_VALID_PREFIX_NAME = \
    And(_alphanum_check(whitelist="._-*"),
        Not(StringIn(["Options"] + [(s + "Reads") for s in _READ_TYPES])))

# Valid paths for prefixes; avoids some problems with e.g. Bowtie2
_VALID_PREFIX_PATH = \
    And(IsStr(), Not(ValuesIntersect("\\:?\"<>|() \t\n\v\f\r")),
        default=REQUIRED_VALUE)

# Valid strings for targets / samples / libraries / lanes
_VALID_TARGET_NAME = \
    And(_alphanum_check(whitelist="._-"),
        ValueGE(2, key=len, description="at least two characters long"))


_VALIDATION_OPTIONS = {
    # Sequencing platform, used to tag read-groups.
    "Platform": StringIn(("CAPILLARY", "LS454", "ILLUMINA", "SOLID",
                          "HELICOS", "IONTORRENT", "PACBIO"),
                         default="ILLUMINA"),
    # Offset for quality scores in FASTQ files.
    "QualityOffset": ValueIn((33, 64, "Solexa"),
                             default=33),
    # Split a lane into multiple entries, one for each (pair of) file(s)
    "SplitLanesByFilenames": Or(IsBoolean, IsListOf(IsStr),
                                default=True),
    # Format to use when compressing FASTQ files ("gz" or "bz2")
    "CompressionFormat": ValueIn(("gz", "bz2"),
                                 default="bz2"),

    "AdapterRemoval": {
        "Version": ValueIn(("v1.4", "v1.5+"),
                           default="v1.5+"),
        "--pcr1": IsStr,
        "--pcr2": IsStr,
        "--adapter1": IsStr,
        "--adapter2": IsStr,
        "--maxns": IsUnsignedInt,
        "--minquality": IsUnsignedInt,
        "--trimns": Or(IsNone, IsBoolean),
        "--trimqualities": Or(IsNone, IsBoolean),
        "--collapse": Or(IsNone, IsBoolean, default=True),
        "--mm": Or(IsFloat, IsUnsignedInt,
                   default=3),
        "--minlength": IsUnsignedInt(default=25),
        "--minalignmentlength": IsUnsignedInt,
        "--shift": IsUnsignedInt,
        "--5prime": IsStr,
        },

    # Which aliger/mapper to use (BWA/Bowtie2)
    "Aligners": {
        "Program": ValueIn(("BWA", "Bowtie2"),
                           default="BWA"),
        "BWA": {
            # Mapping algorithm; availability depends on BWA version
            "Algorithm": StringIn(("backtrack", "mem", "bwasw"),
                                  default="backtrack"),

            # Minimum mapping quality (PHREAD) of reads to retain
            "MinQuality": IsUnsignedInt(default=0),
            # Remove unmapped reads or not
            "FilterUnmappedReads": IsBoolean(default=True),
            # Use seed region during mapping
            # Verbose name for command-line option "-l 65535"
            "UseSeed": IsBoolean(default=True),
            # Any number of user specific options
            StringStartsWith("-"): Or(IsListOf(IsStr, IsInt, IsFloat),
                                      Or(IsStr, IsInt, IsFloat, IsNone)),
        },
        "Bowtie2": {
            # Minimum mapping quality (PHREAD) of reads to retain
            "MinQuality": IsUnsignedInt(default=0),
            # Remove unmapped reads or not
            "FilterUnmappedReads": IsBoolean(default=True),
            # Any number of user specific options
            StringStartsWith("-"): Or(IsListOf(IsStr, IsInt, IsFloat),
                                      Or(IsStr, IsInt, IsFloat, IsNone)),
        },
    },

    # Does sample contain PCR duplicates / what to do about it.
    # True is equivalent of 'remove'.
    "PCRDuplicates": StringIn((True, False, 'mark', 'filter'),
                              default='filter'),
    # Qualities should be rescaled using mapDamage
    "RescaleQualities": IsBoolean(default=False),

    "mapDamage": {
        # Tabulation options
        "--downsample": Or(IsUnsignedInt, IsFloat),
        "--length": IsUnsignedInt,
        "--around": IsUnsignedInt,
        "--min-basequal": IsUnsignedInt,

        # Plotting options
        "--ymax": IsFloat,
        "--readplot": IsUnsignedInt,
        "--refplot": IsUnsignedInt,

        # Model options
        "--rand": IsUnsignedInt,
        "--burn": IsUnsignedInt,
        "--adjust": IsUnsignedInt,
        "--iter": IsUnsignedInt,
        "--forward": IsNone,
        "--reverse": IsNone,
        "--var-disp": IsNone,
        "--jukes-cantor": IsNone,
        "--diff-hangs": IsNone,
        "--fix-nicks": IsNone,
        "--use-raw-nick-freq": IsNone,
        "--single-stranded": IsNone,
        "--seq-length": IsUnsignedInt,
    },

    # Exclude READ_TYPES from alignment/analysis
    "ExcludeReads": Or(IsNone,
                       ValuesSubsetOf(_READ_TYPES),
                       default=[]),

    # Features of pipeline
    "Features": Or(IsNone,
                   ValuesSubsetOf(("Raw BAM", "Realigned BAM", "Coverage",
                                   "Summary", "mapDamage", "Depths",
                                   "DuplicateHist")),
                   default=["Realigned BAM", "Coverage",
                            "Summary", "mapDamage", "Depths"]),
}


_VALIDATION = {
    "Options": _VALIDATION_OPTIONS,

    "Prefixes": {
        _VALID_PREFIX_NAME: {
            "Path": _VALID_PREFIX_PATH,
            "Label": ValueIn(("nuclear", "mitochondrial", "chloroplast",
                              "plasmid", "bacterial", "viral")),
            "RegionsOfInterest": IsDictOf(IsStr, IsStr),
        },
    },

    _VALID_TARGET_NAME: {  # Target
        _VALID_TARGET_NAME: {  # Sample
            _VALID_TARGET_NAME: {  # Library
                _VALID_TARGET_NAME: Or(IsStr, IsDictOf(IsStr, IsStr)),

                "Options": WithoutDefaults(_VALIDATION_OPTIONS),
            },

            "Options": WithoutDefaults(_VALIDATION_OPTIONS),
        },

        "Options": WithoutDefaults(_VALIDATION_OPTIONS),
    },
}


def _mangle_makefile(makefile):
    makefile = copy.deepcopy(makefile)
    makefile["Options"] = makefile["Makefile"].pop("Options")
    makefile["Prefixes"] = makefile["Makefile"].pop("Prefixes")
    makefile["Targets"] = makefile.pop("Makefile")

    _update_options(makefile)
    _update_prefixes(makefile)
    _update_lanes(makefile)
    _update_tags(makefile)

    _split_lanes_by_filenames(makefile)

    return makefile


def _update_options(makefile):
    def _update_possibly_empty_lists(options):
        for key in ("Features", "ExcludeReads"):
            if options[key] is None:
                options[key] = []

    def _do_update_options(options, data, path):
        options = copy.deepcopy(options)
        if "Options" in data:
            if "Features" in data["Options"]:
                raise MakefileError("Features may only be specified at root "
                                    "level, not at %r" % (":".join(path),))

            # Fill out missing values using those of prior levels
            options = fill_dict(destination=data.pop("Options"),
                                source=options)
            _update_possibly_empty_lists(options)

        if len(path) < 2:
            for key in data:
                if key != "Options":
                    _do_update_options(options, data[key], path + (key,))
        else:
            data["Options"] = options

    for data in makefile["Targets"].itervalues():
        _update_possibly_empty_lists(makefile["Options"])
        _do_update_options(makefile["Options"], data, ())


def _update_prefixes(makefile):
    prefixes = {}
    for (name, values) in makefile.get("Prefixes", {}).iteritems():
        filename = values["Path"]
        if "*" in name[:-1]:
            raise MakefileError("The character '*' is not allowed in Prefix "
                                "names; if you use to select .fasta files "
                                "using a search-string, then use the prefix "
                                "name '%s*' instead and specify the wildcards "
                                "in the 'Path' instead."
                                % (name.replace("*", "",)))
        elif name.endswith("*"):
            records = []
            for fname in glob.glob(filename):
                name = os.path.basename(fname).split(".")[0]
                _VALID_PREFIX_NAME(("Prefixes", name), name)
                new_prefix = copy.copy(values)
                new_prefix["Path"] = fname

                records.append((name, new_prefix))
            if not records:
                raise MakefileError("Did not find any matches for glob %s"
                                    % repr(filename))
        else:
            records = [(name, values)]

        for (name, record) in records:
            if name in prefixes:
                raise MakefileError("Multiple prefixes with the same name: %s"
                                    % name)

            if not record["Path"].endswith(".fasta"):
                raise MakefileError("Path for prefix %r does not end with "
                                    ".fasta:\n   %r" % (name, record["Path"]))

            record["Name"] = name
            record["Reference"] = record["Path"]
            prefixes[name] = record

    if not prefixes:
        raise MakefileError("At least one prefix must be specified")
    makefile["Prefixes"] = prefixes


def _update_lanes(makefile):
    prefixes = makefile["Prefixes"]
    for (target_name, samples) in makefile["Targets"].iteritems():
        for (sample_name, libraries) in samples.iteritems():
            for (library_name, lanes) in libraries.iteritems():
                options = lanes.pop("Options")

                for (lane, data) in lanes.iteritems():
                    path = (target_name, sample_name, library_name, lane)
                    lane_type = _determine_lane_type(prefixes, data, path)
                    lanes[lane] = {"Type": lane_type,
                                   "Data": data,
                                   "Options": options}


def _determine_lane_type(prefixes, data, path):
    if isinstance(data, types.StringTypes):
        return "Raw"
    elif isinstance(data, types.DictType):
        if all((key in _READ_TYPES) for key in data):
            for (key, files) in data.iteritems():
                is_paired = paths.is_paired_end(files)

                if is_paired and (key != "Paired"):
                    raise MakefileError("Error at Barcode level; Path "
                                        "includes {Pair} key, but read-type "
                                        "is not Paired:\n    "
                                        "%s:%s" % (":".join(path), key))
                elif not is_paired and (key == "Paired"):
                    raise MakefileError("Error at Barcode level; Paired pre-"
                                        "trimmed reads specified, but path "
                                        "does not contain {Pair} key:\n    "
                                        "%s:%s" % (":".join(path), key))

            return "Trimmed"
        elif all((key in prefixes) for key in data):
            return "BAMs"

    raise MakefileError("Error at Barcode level; keys must either be "
                        "prefix-names, OR 'Paired', 'Single' or 'Collapsed'. "
                        "Found: %s" % (", ".join(data),))


def _update_tags(makefile):
    for (target, samples) in makefile["Targets"].iteritems():
        for (sample, libraries) in samples.iteritems():
            for (library, barcodes) in libraries.iteritems():
                for (barcode, record) in barcodes.iteritems():
                    tags = {"Target": target,
                            "ID": library,
                            "SM": sample,
                            "LB": library,
                            # Source/Current PU may differ if a lane has been
                            # split by filenames, in which case PU_src contains
                            # the original PU, and PU_cur is a derived PU.
                            "PU_src": barcode,
                            "PU_cur": barcode,
                            "PG": record["Options"]["Aligners"]["Program"],
                            "PL": record["Options"]["Platform"].upper()}

                    record["Tags"] = tags


def _split_lanes_by_filenames(makefile):
    iterator = _iterate_over_records(makefile)
    for (target, sample, library, barcode, record) in iterator:
        if record["Type"] == "Raw":
            template = record["Data"]
            record["Data"] = files = paths.collect_files(template)
            split = record["Options"]["SplitLanesByFilenames"]

            if (split == True) or (isinstance(split, list) and (barcode in split)):
                if any(missing_files(file_set) for file_set in files.itervalues()):
                    raise MakefileError("Unable to split by filename for "
                                        "search-string '%s', did not find any "
                                        "files; please verify that the path"
                                        "is correct and update the makefile."
                                        % template)
                elif any(len(v) > 1 for v in files.itervalues()):
                    template = makefile["Targets"][target][sample][library].pop(barcode)
                    keys = ("SE",) if ("SE" in files) else ("PE_1", "PE_2")

                    input_files = [files[key] for key in keys]
                    input_files_iter = itertools.izip_longest(*input_files)
                    for (index, filenames) in enumerate(input_files_iter, start=1):
                        assert len(filenames) == len(keys)
                        assert len(filenames[0]) == len(filenames[-1])
                        new_barcode = "%s_%03i" % (barcode, index)

                        current = copy.deepcopy(template)
                        current["Data"] = dict((key, [filename]) for (key, filename) in zip(keys, filenames))
                        current["Tags"]["PU_cur"] = new_barcode

                        makefile["Targets"][target][sample][library][new_barcode] = current


def _validate_makefiles(config, makefiles):
    for makefile in makefiles:
        _validate_makefile_libraries(makefile)
        _validate_makefile_adapters(makefile)
    _validate_makefiles_duplicate_targets(config, makefiles)
    _validate_makefiles_duplicate_files(makefiles)
    _validate_makefiles_features(makefiles)
    _validate_hg_prefixes(makefiles)

    return makefiles


def _validate_makefile_adapters(makefile):
    """Checks for the default adapter sequences specified in the wrong
    orientation for AdapterRemoval, which is a typical mistake when using
    the --pcr2 option.
    """
    # The non-reverse complemented mate 2 adapter, as seen in raw FASTQ reads
    adapter_2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT"

    tests = {
        # --pcr2 expects the reverse complement of the mate 2 adapter seq.
        "--pcr2": adapter_2,
        # --adapter2 (AdapterRemoval v2) expects the regular sequence
        "--adapter2": sequences.reverse_complement(adapter_2)
    }

    def check_options(options, results):
        for key, value in tests.iteritems():
            if options.get(key) == value:
                results[key] = True

    results = dict.fromkeys(tests, False)
    for (_, _, _, _, record) in _iterate_over_records(makefile):
        adapterrm_opt = record.get("Options", {}).get("AdapterRemoval", {})
        check_options(adapterrm_opt, results)

    adapterrm_opt = makefile.get("Options", {}).get("AdapterRemoval", {})
    check_options(adapterrm_opt, results)

    if any(results.itervalues()):
        print_warn("WARNING: An adapter specified for AdapterRemoval "
                   "corresponds to the default sequence, but is reverse "
                   "complemented. Please make sure that this is intended! ",
                   end="")

        if results["--pcr2"]:
            print_warn("For --pcr2, the sequence given should be the "
                       "reverse complement of the sequence observed in the "
                       "mate 2 FASTQ file.\n")

        if results["--adapter2"]:
            print_warn("For --adapter2 (AdapterRemoval v2, only) the value "
                       "should be exactly as observed in the FASTQ reads.\n")


def _validate_makefile_libraries(makefile):
    libraries = collections.defaultdict(set)
    iterator = _iterate_over_records(makefile)
    for (target, sample, library, _, _) in iterator:
        libraries[(target, library)].add(sample)

    for ((target, library), samples) in libraries.iteritems():
        if len(samples) > 1:
            raise MakefileError("Library '%s' in target '%s' spans multiple "
                                " samples: %s" % (library, target,
                                                  ", ".join(samples)))


def _validate_makefiles_duplicate_files(makefiles):
    filenames = collections.defaultdict(list)
    for makefile in makefiles:
        iterator = _iterate_over_records(makefile)
        for (target, sample, library, barcode, record) in iterator:
            current_filenames = []
            if record["Type"] == "Raw":
                for raw_filenames in record["Data"].itervalues():
                    current_filenames.extend(raw_filenames)
            else:
                current_filenames.extend(record["Data"].values())

            for realpath in map(os.path.realpath, current_filenames):
                filenames[realpath].append((target, sample, library, barcode))

    has_overlap = {}
    for (filename, records) in filenames.iteritems():
        if len(records) > 1:
            has_overlap[filename] = list(set(records))

    by_records = sorted(zip(has_overlap.values(), has_overlap.keys()))
    for (records, pairs) in itertools.groupby(by_records, lambda x: x[0]):
        pairs = list(pairs)
        description = _describe_files_in_multiple_records(records, pairs)

        if len(set(record[0] for record in records)) != len(records):
            message = "Path included multiple times in target:\n"
            raise MakefileError(message + description)
        else:
            print_warn("WARNING: Path included in multiple targets:",
                       file=sys.stderr)
            print_warn(description, file=sys.stderr)
            print_warn(file=sys.stderr)


def _describe_files_in_multiple_records(records, pairs):
    descriptions = []
    for (index, record) in enumerate(records, start=1):
        descriptions.append("\t- Record {0}: Name: {1},  Sample: {2},  "
                            "Library: {3},  Barcode: {4}".format(index,
                                                                 *record))

    for (index, (_, filename)) in enumerate(sorted(pairs), start=1):
        message = "\t- Canonical path {0}: {1}"
        descriptions.append(message.format(index, filename))

    return "\n".join(descriptions)


def _validate_makefiles_duplicate_targets(config, makefiles):
    targets = set()
    for makefile in makefiles:
        destination = config.destination
        if destination is None:
            filename = makefile["Statistics"]["Filename"]
            destination = os.path.dirname(filename)

        for target in makefile["Targets"]:
            key = (destination, target)
            if key in targets:
                raise MakefileError("Target name '%s' used multiple times; "
                                    "output files would be clobbered!"
                                    % target)
            targets.add(key)


def _validate_makefiles_features(makefiles):
    for makefile in makefiles:
        features = makefile["Options"]["Features"]
        roi_enabled = False

        for prefix in makefile["Prefixes"].itervalues():
            roi_enabled |= bool(prefix.get("RegionsOfInterest"))

        if "Depths" in features and roi_enabled:
            if not (("Raw BAM" in features) or ("Realigned BAM") in features):
                raise MakefileError("The feature 'Depths' (depth histograms) "
                                    "with RegionsOfInterest enabled, requires "
                                    "that either the feature 'Raw BAM' or the "
                                    "feature 'Raligned BAM' is enabled.")


def _validate_hg_prefixes(makefiles):
    """Implementation of the checks included in GATK, which require that the
    FASTA for the human genome is ordered 1 .. 23, .
    """
    already_validated = set()
    print_info("  - Validating prefixes ...", file=sys.stderr)
    for makefile in makefiles:
        uses_gatk = "Realigned BAM" in makefile["Options"]["Features"]
        for prefix in makefile["Prefixes"].itervalues():
            path = prefix["Path"]
            if path in already_validated:
                continue

            if not os.path.exists(path):
                print_info("    - Reference FASTA file does not exist:\n"
                           "      %r" % (path,), file=sys.stderr)
                continue
            elif not os.path.exists(path + ".fai"):
                print_info("    - Index does not exist for %r; this may "
                           "take a while ..." % (path,), file=sys.stderr)

                if not os.access(os.path.dirname(path), os.W_OK):
                    message = \
                        "FASTA index for prefix is missing, but folder is\n" \
                        "not writable, so it cannot be created:\n" \
                        "  Prefix = %s\n\n" \
                        "Either change permissions on the folder, or move\n" \
                        "the prefix to different location." % (path,)
                    raise MakefileError(message)

                # Use pysam to index the file
                pysam.Fastafile(path).close()

            contigs = []
            with open(path + ".fai") as handle:
                for line in handle:
                    name, size, _ = line.split('\t', 2)
                    contigs.append((name, int(size)))

            _do_validate_hg_prefix(makefile, prefix, contigs, fatal=uses_gatk)
            already_validated.add(path)


def _do_validate_hg_prefix(makefile, prefix, contigs, fatal):
    if not _is_invalid_hg_prefix(contigs):
        return

    message = \
        "Prefix appears to be a human genome, but chromosomes are ordered\n" \
        "lexically (chr1, chr10, chr11, ...), rather than numerically\n" \
        "(chr1, chr2, chr3, ...):\n\n" \
        "  Makefile = %s\n" \
        "  Prefix   = %s\n\n" \
        "GATK requires that human chromosomes are ordered numerically;\n%s\n" \
        "See the documentation at the GATK website for more information:\n  " \
        "http://www.broadinstitute.org/gatk/guide/article?id=1204\n"

    prefix_path = prefix["Path"]
    mkfile_path = makefile["Statistics"]["Filename"]
    if fatal:
        details = "Either disable GATK in the makefile, or fix the prefix."
        message %= (mkfile_path, prefix_path, details)

        raise MakefileError(message)
    else:
        details = \
            "You will not be able to use the resulting BAM file with GATK."
        message %= (mkfile_path, prefix_path, details)
        print_warn("\nWARNING:\n", message, file=sys.stderr, sep="")


def _is_invalid_hg_prefix(contigs):
    hg_contigs = {
        # Contig sizes based on hg18 and hg19 and hg38
        "chr1": [247249719, 249250621, 248956422],
        "chr2": [242951149, 243199373, 242193529],
        "chr10": [135374737, 135534747, 133797422],
    }

    size_to_idx = dict((size, idx) for (idx, (_, size)) in enumerate(contigs))

    # Equivalent to the GATK 'nonCanonicalHumanContigOrder' function
    for (key, values) in hg_contigs.iteritems():
        for value in values:
            if value in size_to_idx:
                hg_contigs[key] = size_to_idx[value]
                break
        else:
            # Contig not found; probably not hg18, hg19, or hg38
            return False

    return not (hg_contigs["chr1"] < hg_contigs["chr2"] < hg_contigs["chr10"])


def _iterate_over_records(makefile):
    for (target, samples) in makefile["Targets"].items():
        for (sample, libraries) in samples.items():
            for (library, barcodes) in libraries.items():
                for (barcode, record) in barcodes.items():
                    yield target, sample, library, barcode, record
