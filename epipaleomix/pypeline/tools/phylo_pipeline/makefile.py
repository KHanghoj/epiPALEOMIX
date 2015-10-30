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
import types
import pysam

import pypeline.common.makefile
from pypeline.common.makefile import \
    MakefileError, \
    REQUIRED_VALUE, \
    IsDictOf, \
    IsListOf, \
    IsInt, \
    IsStr, \
    StringIn, \
    IsFloat, \
    IsUnsignedInt, \
    IsBoolean, \
    IsNone, \
    ValueIn, \
    ValuesSubsetOf, \
    StringStartsWith, \
    StringEndsWith, \
    CLI_PARAMETERS, \
    And, \
    Or, \
    Not

from pypeline.common.fileutils import \
    swap_ext, \
    add_postfix
from pypeline.common.utilities import \
    fill_dict
from pypeline.common.console import \
    print_info, \
    print_warn
from pypeline.common.text import \
    parse_padded_table
from pypeline.common.bedtools import \
    BEDRecord


def read_makefiles(options, filenames, commands):
    print_info("Reading makefile(s):")
    steps = frozenset(key for (key, _) in commands)

    makefiles = []
    for filename in filenames:
        makefile = pypeline.common.makefile.read_makefile(filename, _VALIDATION)
        makefile = _mangle_makefile(options, makefile["Makefile"], steps)
        makefiles.append(makefile)
    return makefiles


def _mangle_makefile(options, mkfile, steps):
    _collapse_samples(mkfile)
    _update_regions(options, mkfile)
    _update_subsets(mkfile, steps)
    _update_filtering(mkfile)
    _update_sample_sets(mkfile)
    _update_genotyping(mkfile)
    _update_msa(mkfile)
    _update_homozygous_contigs(mkfile)
    _check_bam_sequences(options, mkfile, steps)
    _check_genders(mkfile)
    _update_and_check_max_read_depth(options, mkfile)
    _check_indels_and_msa(mkfile)
    mkfile["Nodes"] = ()

    return mkfile


def _collapse_samples(mkfile):
    groups, samples = {}, set()
    def _collect_samples(samples_dict, path = ()):
        current_samples = {}
        for (key, subdd) in samples_dict.iteritems():
            if key.startswith("<") and key.endswith(">"):
                key = key.lstrip("<").rstrip(">")
                current_samples.update(_collect_samples(subdd, path + (key,)))
            elif key not in samples:
                samples.add(key)
                subdd["Name"] = key
                current_samples[key] = subdd
            else:
                raise MakefileError("Duplicate sample-name: %r" % (key,))

        groups[path] = current_samples
        return current_samples

    _collect_samples(mkfile["Project"]["Samples"])
    mkfile["Project"]["Samples"] = groups.pop(())
    mkfile["Project"]["Groups"] = groups


def _select_samples(select, groups, samples, path):
    selection = set()
    for group in select:
        if group.startswith("<") and group.endswith(">"):
            key = tuple(group[1:-1].split("/"))
            if key not in groups:
                raise MakefileError("Unknown group specifed for filtering %r: %r" % (path, key))
            selection.update(groups[key])
        elif group in samples:
            selection.add(group)
        else:
            raise MakefileError("Unknown/Invalid group specifed for filtering %r: %r" % (path, group))
    return selection


def _update_regions(options, mkfile):
    print_info("    - Validating regions of interest ...")
    mkfile["Project"]["Regions"] = mkfile["Project"].pop("RegionsOfInterest")

    for (name, subdd) in mkfile["Project"]["Regions"].iteritems():
        if "Prefix" not in subdd:
            raise MakefileError("No genome specified for regions %r" % (name,))

        subdd["Name"]   = name
        subdd["Desc"]   = "{Prefix}.{Name}".format(**subdd)
        subdd["BED"]    = os.path.join(options.regions_root, subdd["Desc"] + ".bed")
        subdd["FASTA"]  = os.path.join(options.prefix_root, subdd["Prefix"] + ".fasta")

        required_files = (
            ("Regions file", subdd["BED"], None),
            ("Reference sequence", subdd["FASTA"], None),
            ("Reference sequence index", subdd["FASTA"] + ".fai",
             "Please index using 'samtools faidx %s'" % (subdd["FASTA"],)))

        for (desc, path, instructions) in required_files:
            if not os.path.isfile(path):
                message = "%s does not exist for %r:\n  Path = %r" \
                                % (desc, name, path)
                if instructions:
                    message = "%s\n%s" % (message, instructions)
                raise MakefileError(message)

        # Collects seq. names / validate regions
        subdd["Sequences"] = {None : _collect_and_validate_regions(subdd)}
        subdd["SubsetFiles"] = {None : ()}

        sampledd = subdd["Genotypes"] = {}
        for sample_name in mkfile["Project"]["Samples"]:
            fasta_file = ".".join((sample_name, subdd["Desc"], "fasta"))
            sampledd[sample_name] = os.path.join(options.destination,
                                              mkfile["Project"]["Title"],
                                              "genotypes",
                                              fasta_file)


_CONTIGS_CACHE = {}
def _collect_fasta_contigs(regions):
    filename = regions["FASTA"] + ".fai"
    if filename in _CONTIGS_CACHE:
        return _CONTIGS_CACHE[filename]

    contigs = {}
    with open(filename) as faihandle:
        for line in faihandle:
            name, length, _ = line.split(None, 2)
            if name in contigs:
                message = ("Reference contains multiple identically named "
                           "sequences:\n  Path = %r\n  Name = %r\n"
                           "Please ensure that sequences have unique names") \
                           % (regions["FASTA"], name)
                raise MakefileError(message)

            contigs[name] = int(length)

    _CONTIGS_CACHE[filename] = contigs
    return contigs


def _collect_and_validate_regions(regions):
    contigs = _collect_fasta_contigs(regions)
    sequences = set()
    with open(regions["BED"]) as bedhandle:
        for (line_num, line) in enumerate(bedhandle):
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            try:
                bed = BEDRecord(line)
            except ValueError, error:
                raise MakefileError(("Error parsing line %i in regions file:\n"
                                     "  Path = %r\n  Line = %r\n\n%s")
                                    % (line_num + 1, regions["BED"],
                                       line, error))

            if len(bed) < 6:
                url = "http://genome.ucsc.edu/FAQ/FAQformat.html#format1"
                name = repr(bed.name) if len(bed) > 3 else "unnamed record"
                raise MakefileError(("Region at line #%i (%s) does not "
                                     "contain the expected number of fields; "
                                     "the first 6 fields are required. C.f. "
                                     "defination at\n   %s\n\nPath = %r")
                                    % (line_num, name, url, regions["BED"]))

            contig_len = contigs.get(bed.contig)
            if contig_len is None:
                raise MakefileError(("Regions file contains contig not found "
                                     "in reference:\n  Path = %r\n  Contig = "
                                     "%r\n\nPlease ensure that all contig "
                                     "names match the reference names!")
                                    % (regions["BED"], bed.contig))
            elif not (0 <= bed.start < bed.end <= contig_len):
                raise MakefileError(("Regions file contains invalid region:\n"
                                     "  Path   = %r\n  Contig = %r\n"
                                     "  Start  = %s\n  End    = %s\n\n"
                                     "Expected 0 <= Start < End <= %i!")
                                    % (regions["BED"], bed.contig, bed.start,
                                       bed.end, contig_len))

            sequences.add(bed.name)

    return frozenset(sequences)


def _update_subsets(mkfile, steps):
    subsets_by_regions = mkfile["Project"]["Regions"]
    def _collect_subsets(roi, subset, path):
        if roi not in subsets_by_regions:
            raise MakefileError("Subset of unknown region (%r) requested at %r" % (roi, path))

        roi_fname = swap_ext(subsets_by_regions[roi]["BED"], subset + ".names")
        if not os.path.isfile(roi_fname):
            raise MakefileError(("Subset file does not exist for Regions Of Interest:\n"
                                 "  Region = %r\n  Subset = %r\n  Path   = %r")
                                 % (roi, subset, roi_fname))

        sequences = set()
        with open(roi_fname) as handle:
            for line in handle:
                line = line.strip()
                if line and not line.startswith("#"):
                    sequences.add(line)

        known_seqs   = subsets_by_regions[roi]["Sequences"][None]
        unknown_seqs = sequences - known_seqs
        if unknown_seqs:
            message = ("Unknown sequences in subset file:\n"
                       "  File   = %r\n  Region = %r\n  Subset = %r\n"
                       "  Unknown sequence names =") \
                       % (roi_fname, roi, subset)
            unknown_seqs = list(sorted(unknown_seqs))
            if len(unknown_seqs) > 5:
                unknown_seqs = unknown_seqs[:5]  + ["..."]
            message = "\n    - ".join([message] + unknown_seqs)
            raise MakefileError(message)

        subsets_by_regions[roi]["SubsetFiles"][subset] = (roi_fname,)
        subsets_by_regions[roi]["Sequences"][subset] = frozenset(sequences)

    if "phylogeny:examl" in steps:
        for (key, subdd) in mkfile["PhylogeneticInference"].iteritems():
            for (subkey, roidd) in subdd["RegionsOfInterest"].iteritems():
                if subkey not in subsets_by_regions:
                    message = \
                        "Unknown regions name in phylogenetic inference:\n" \
                        "\tPath = PhylogeneticInference:%s:RegionsOfInterest" \
                        "\n\tName = %s"
                    raise MakefileError(message % (key, subkey))

                roidd["Name"] = subkey

                if roidd.get("SubsetRegions") is not None:
                    path = "PhylogeneticInference:%s:RegionsOfInterest:%s" % (key, subkey)
                    _collect_subsets(subkey, roidd["SubsetRegions"], path)

    if "paml:codeml" in steps:
        for (roi, subset) in mkfile["PAML"]["codeml"]["SubsetRegions"].iteritems():
            _collect_subsets(roi, subset, "PAML:codeml:SubsetRegions")


def _update_filtering(mkfile):
    samples = mkfile["Project"]["Samples"]
    groups  = mkfile["Project"]["Groups"]

    filtering = {}
    for (target, filter_by) in mkfile["Project"]["FilterSingletons"].iteritems():
        if target.startswith("<") and target.endswith(">"):
            raise MakefileError("Singleton-filtering must be specified per "
                                "sample, not by groups: %r" % (target,))
        elif target not in samples:
            raise MakefileError("Unknown/Invalid sample specifed for singleton filtering: %r" % (target,))
        elif target in filter_by:
            raise MakefileError("Attempting to filter singleton in sample using itself as comparison: %r" % (target,))

        path = "Project:FilterSingletons:%s" % (target,)
        filtering[target] = _select_samples(filter_by, groups, samples, path)

        # Implicit inclusion is allowed, since that is useful in some cases,
        # where we want to filter a sample based on the group it is a member of
        if target in filtering[target]:
            # The target itself must be excluded, as including it is invalid
            filtering[target] = filtering[target] - set((target,))
            print_warn("Warning: Sample %r is singleton-filtered using a "
                       "group it is also a member of; this may be by mistake."
                       % (target,))

        if not filtering[target]:
            raise MakefileError("No samples specified by which to "
                                "singleton-filter by for %r" % (target,))

    mkfile["Project"]["FilterSingletons"] = filtering


def _update_homozygous_contigs(mkfile):
    """Treat unspecified values for HomozygousContigs as an empty list, in
    order that the user does not need to specify "[]" for empty lists.
    """
    for regions in mkfile["Project"]["Regions"].itervalues():
        hcontigs = regions["HomozygousContigs"]
        for key, contigs in hcontigs.items():
            if contigs is None:
                hcontigs[key] = []


def _check_bam_sequences(options, mkfile, steps):
    """Check that the BAM files contains the reference sequences found in the
    FASTA file, matched by name and length; extra sequences are permitted. This
    check is only done if genotyping is to be carried out, to reduce the
    overhead of reading the BAM file headers.

    """
    if ("genotype" not in steps) and ("genotyping" not in steps):
        return

    print_info("    - Validating BAM files ...")
    bam_files = {}
    for regions in mkfile["Project"]["Regions"].itervalues():
        for sample in mkfile["Project"]["Samples"].itervalues():
            filename = os.path.join(options.samples_root, "%s.%s.bam"
                                    % (sample["Name"], regions["Prefix"]))
            if regions["Realigned"]:
                filename = add_postfix(filename, ".realigned")

            if os.path.exists(filename):
                bam_files[filename] = _collect_fasta_contigs(regions)

    for (filename, contigs) in bam_files.iteritems():
        with pysam.Samfile(filename) as handle:
            bam_contigs = dict(zip(handle.references, handle.lengths))

            for (contig, length) in contigs.iteritems():
                bam_length = bam_contigs.get(contig)

                if bam_length is None:
                    message = ("Reference sequence missing from BAM file; "
                               "BAM file aligned against different prefix?\n"
                               "    BAM file = %s\n    Sequence name = %s") \
                               % (filename, contig)
                    raise MakefileError(message)
                elif bam_length != length:
                    message = ("Length of reference sequence in FASTA differs "
                               "from length of sequence in BAM file; BAM file "
                               "aligned against different prefix?\n"
                               "    BAM file = %s\n"
                               "    Length in FASTA = %s\n"
                               "    Length in BAM = %s") \
                               % (filename, length, bam_length)
                    raise MakefileError(message)


def _check_genders(mkfile):
    all_contigs = set()
    contigs_genders = set()
    regions_genders = set()
    for regions in mkfile["Project"]["Regions"].itervalues():
        all_contigs.update(_collect_fasta_contigs(regions))

        for contigs in regions["HomozygousContigs"].itervalues():
            contigs_genders.update(contigs)

        current_genders = set(regions["HomozygousContigs"])
        if not regions_genders:
            regions_genders = current_genders
        elif regions_genders != current_genders:
            raise MakefileError("List of genders for regions %r does not "
                                "match other regions" % (regions["Name"],))

    if not regions_genders:
        raise MakefileError("No genders have been specified in makefile; "
                            "please list all sample genders and assosiated "
                            "homozygous contigs (if any).")

    for sample in mkfile["Project"]["Samples"].itervalues():
        if sample["Gender"] not in regions_genders:
            genders = ", ".join(map(repr, regions_genders))
            message = "Sample %r has unknown gender %r; known genders are %s" \
                % (sample["Name"], sample["Gender"], genders)
            raise MakefileError(message)

    unknown_contigs = contigs_genders - all_contigs
    if unknown_contigs:
        print_warn("WARNING: Unknown contig(s) in 'HomozygousContigs':\n    - "
                   + "\n    - ".join(unknown_contigs))
        print_warn("Please verify that the list(s) of contigs is correct!")


def _update_and_check_max_read_depth(options, mkfile):
    if any(subdd["VCF_Filter"]["MaxReadDepth"] == "auto"
           for subdd in mkfile["Genotyping"].itervalues()):
        print_info("    - Determinining max-depth from depth-histograms ...")

    for (key, settings) in mkfile["Genotyping"].iteritems():
        required_keys = set()
        for sample in mkfile["Project"]["Samples"].itervalues():
            if sample["GenotypingMethod"].lower() == "samtools":
                required_keys.add(sample["Name"])

        max_depths = settings["VCF_Filter"]["MaxReadDepth"]
        if isinstance(max_depths, types.DictType):
            # Extra keys are allowed, to make it easier
            # to temporarily disable a sample
            missing_keys = required_keys - set(max_depths)
            if missing_keys:
                missing_keys = "\n    - ".join(sorted(missing_keys))
                message = "MaxReadDepth not specified for the following " \
                          "samples for %r:\n    - %s" % (key, missing_keys)
                raise MakefileError(message)

        elif isinstance(max_depths, types.StringTypes):
            assert max_depths.lower() == "auto", max_depths
            prefix = mkfile["Project"]["Regions"][key]["Prefix"]
            max_depths = {}

            for sample in required_keys:
                fname = "%s.%s.depths" % (sample, prefix)
                fpath = os.path.join(options.samples_root, fname)
                max_depths[sample] = _read_max_depths(fpath, prefix, sample)

            settings["VCF_Filter"]["MaxReadDepth"] = max_depths
        else:
            max_depths = dict.fromkeys(required_keys, max_depths)
            settings["VCF_Filter"]["MaxReadDepth"] = max_depths


def _read_max_depths(filename, prefix, sample):
    if filename in _DEPTHS_CACHE:
        return _DEPTHS_CACHE[filename]

    max_depth = None
    try:
        with open(filename) as handle:
            for row in parse_padded_table(handle):
                if row["Name"] == sample and \
                        row["Sample"] == "*" and \
                        row["Library"] == "*" and \
                        row["Contig"] == "*":
                    max_depth = row["MaxDepth"]
                    break
            else:
                raise MakefileError("Could not find MaxDepth in "
                                    "depth-histogram: %r" % (filename,))

    except (OSError, IOError), error:
        raise MakefileError("Error reading depth-histogram (%s): %s"
                            % (filename, error))

    if max_depth == "NA":
        raise MakefileError("MaxDepth is not calculated for sample (%s);\n"
                            "cannot determine MaxDepth values automatically."
                            % (filename,))
    max_depth = int(max_depth)

    print_info("        - %s.%s = %i" % (sample, prefix, max_depth))
    _DEPTHS_CACHE[filename] = max_depth
    return max_depth


_DEPTHS_CACHE = {}


def _check_indels_and_msa(mkfile):
    msa     = mkfile["MultipleSequenceAlignment"]
    regions = mkfile["Project"]["Regions"]
    for (name, subdd) in regions.iteritems():
        msa_enabled = msa[name]["Enabled"]

        if subdd["IncludeIndels"] and not msa_enabled:
            raise MakefileError("Regions %r includes indels, but MSA is disabled!" % (name,))


def _update_sample_sets(mkfile):
    samples = mkfile["Project"]["Samples"]
    groups  = mkfile["Project"]["Groups"]

    for (key, subdd) in mkfile["PhylogeneticInference"].iteritems():
        subdd["ExcludeSamples"] = \
          _select_samples(subdd["ExcludeSamples"], groups, samples, "PhylogeneticInference:%s:ExcludeSamples" % (key,))

        # Replace None with an empty list, to simplify code using this value
        root_trees_on = subdd["RootTreesOn"] or ()
        subdd["RootTreesOn"] = \
          _select_samples(root_trees_on, groups, samples, "PhylogeneticInference:%s:RootTreesOn" % (key,))

    mkfile["PAML"]["codeml"]["ExcludeSamples"] = \
      _select_samples(mkfile["PAML"]["codeml"]["ExcludeSamples"], groups, samples, "PAML:codeml:ExcludeSamples")


def _update_genotyping(mkfile):
    genotyping = mkfile["Genotyping"]
    defaults   = genotyping.pop("Defaults")
    defaults.setdefault("Padding", 5)
    defaults["VCF_Filter"].setdefault("MaxReadDepth", 0)

    for (key, subdd) in genotyping.iteritems():
        if subdd.get("GenotypeEntirePrefix"):
            message = "GenotypeEntirePrefix is only allowed for prefixes " \
                      "using default parameters, but is set for %r" % (key,)
            raise MakefileError(message)

    for key in mkfile["Project"]["Regions"]:
        subdd = fill_dict(genotyping.get(key, {}), defaults)
        subdd["Random"]["--padding"] = subdd["Padding"]
        genotyping[key] = subdd

    regions = set(genotyping)
    unknown_regions = regions - set(mkfile["Project"]["Regions"])
    if unknown_regions:
        raise MakefileError("Unknown Regions of Interest in Genotyping: %s" \
                            % (", ".join(unknown_regions),))


def _update_msa(mkfile):
    msa      = mkfile["MultipleSequenceAlignment"]
    defaults = msa.pop("Defaults")
    defaults.setdefault("Program", "MAFFT")
    defaults["MAFFT"].setdefault("Algorithm", "MAFFT")

    for key in mkfile["Project"]["Regions"]:
        msa[key] = fill_dict(msa.get(key, {}), defaults)

    unknown_regions = set(msa) - set(mkfile["Project"]["Regions"])
    if unknown_regions:
        raise MakefileError("Unknown Regions of Interest in Genotyping: %s" \
                            % (", ".join(unknown_regions),))


# Recursive definition of sample tree
_VALIDATION_SUBSAMPLE_KEY = And(StringStartsWith("<"),
                                StringEndsWith(">"))
_VALIDATION_SAMPLES_KEY = And(IsStr, Not(_VALIDATION_SUBSAMPLE_KEY))
_VALIDATION_SAMPLES = {
    _VALIDATION_SAMPLES_KEY: {
        "GenotypingMethod": StringIn(("reference sequence",
                                      "random sampling",
                                      "samtools"),
                                     default="samtools"),
        "SpeciesName": IsStr,  # Not used; left for backwards compatibility
        "CommonName": IsStr,   # Not used; left for backwards compatibility
        "Gender": IsStr(default=REQUIRED_VALUE),
    }
}
_VALIDATION_SAMPLES[_VALIDATION_SUBSAMPLE_KEY] = _VALIDATION_SAMPLES

# Genotyping settings; note that explicit lists must not be used here, to allow
# proper inheritance of default values. Use IsListOf instead.
_VALIDATION_GENOTYPES = {
    "Padding": IsUnsignedInt,
    "GenotypeEntirePrefix": IsBoolean(default=False),
    "MPileup": {
        StringStartsWith("-"): Or(IsInt, IsStr, IsNone),
    },
    "BCFTools": {
        StringStartsWith("-"): Or(IsInt, IsStr, IsNone),
    },
    "Random": {
        "--min-distance-to-indels": IsUnsignedInt,
    },
    "VCF_Filter": {
        "MaxReadDepth": Or(IsUnsignedInt, IsDictOf(IsStr, IsUnsignedInt),
                           StringIn(("auto",))),

        "--keep-ambigious-genotypes": IsNone,
        "--min-quality": IsUnsignedInt,
        "--min-allele-frequency": IsFloat,
        "--min-mapping-quality": IsUnsignedInt,
        "--min-read-depth": IsUnsignedInt,
        "--max-read-depth": IsUnsignedInt,
        "--min-num-alt-bases": IsUnsignedInt,
        "--min-distance-to-indels": IsUnsignedInt,
        "--min-distance-between-indels": IsUnsignedInt,
        "--min-strand-bias": IsFloat,
        "--min-baseq-bias": IsFloat,
        "--min-mapq-bias": IsFloat,
        "--min-end-distance-bias": IsFloat,
    },
}

_VALIDATION_MSA = {
    "Enabled": IsBoolean(default=True),
    "Program": StringIn(("mafft",)),  # TODO: Add support for other programs

    "MAFFT": {
        "Algorithm": StringIn(("mafft", "auto",
                               "FFT-NS-1", "FFT-NS-2", "FFT-NS-i",
                               "NW-INS-i", "L-INS-i", "E-INS-i", "G-INS-i")),
        StringStartsWith("-"): CLI_PARAMETERS,
    },
}


_VALIDATION = {
    "Project": {
        "Title": IsStr(default="Untitled"),
        "Samples": _VALIDATION_SAMPLES,
        "RegionsOfInterest": {
            IsStr: {
                "Prefix": IsStr(default=REQUIRED_VALUE),
                "Realigned": IsBoolean(default=False),
                "ProteinCoding": IsBoolean(default=False),
                "IncludeIndels": IsBoolean(default=True),
                "HomozygousContigs": {
                    IsStr: Or(IsNone, IsListOf(IsStr))
                    },
                },
            },
        "FilterSingletons": {
            IsStr: [IsStr],
            },
        },
    "Genotyping": {
        "Defaults": _VALIDATION_GENOTYPES,
        IsStr: _VALIDATION_GENOTYPES,
    },
    "MultipleSequenceAlignment": {
        "Defaults": _VALIDATION_MSA,
        IsStr: _VALIDATION_MSA,
        },
    "PhylogeneticInference": {
        IsStr: {
            # Which program to use; TODO: Add support for other programs
            "Program": StringIn(("examl",), default="examl"),
            # Exclude one or more samples from the phylogeny
            "ExcludeSamples": [IsStr],
            # Which samples to root the final trees on / or midpoint rooting
            "RootTreesOn": [IsStr],
            # Create a tree per gene, for each region of interest,
            # or create a supermatrix tree from all regions specified.
            "PerGeneTrees": IsBoolean(default=False),
            # Selection of regions of interest / settings per region
            "RegionsOfInterest": {
                IsStr: {
                    "Partitions": Or(And(IsStr,
                                         ValuesSubsetOf("123456789X")),
                                     ValueIn([False]),
                                     default=REQUIRED_VALUE),
                    "SubsetRegions": Or(IsStr, IsNone, default=None),
                },
            },
            "SubsetRegions": {
                IsStr: IsStr,
            },
            "ExaML": {
                "Bootstraps": IsUnsignedInt(default=100),
                "Replicates": IsUnsignedInt(default=1),
                "Model": StringIn(("GAMMA", "PSR"),
                                  default="gamma"),
            }
        }
    },
    "PAML": {
        "codeml": {
            "ExcludeSamples": [IsStr],
            "SubsetRegions": {
                IsStr: IsStr,
            },
            IsStr: {
                "ControlFile": IsStr(default=REQUIRED_VALUE),
                "TreeFile": IsStr(default=REQUIRED_VALUE),
            },
        },
    },
}
