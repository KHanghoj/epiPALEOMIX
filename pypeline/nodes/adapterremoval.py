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

from pypeline.node import \
    CommandNode
from pypeline.atomiccmd.command import \
    AtomicCmd,\
    CmdError
from pypeline.atomiccmd.sets import \
    ParallelCmds
from pypeline.atomiccmd.builder import \
    AtomicCmdBuilder, \
    use_customizable_cli_parameters, \
    create_customizable_cli_parameters
from pypeline.nodes.validation import \
    check_fastq_files

import pypeline.common.utilities as utilities
import pypeline.common.fileutils as fileutils
import pypeline.common.versions as versions
import pypeline.tools.factory as factory


VERSION_14 = "1.4"
_VERSION_14_CHECK = versions.Requirement(call=("AdapterRemoval", "--version"),
                                         search=r"ver. (\d+)\.(\d+)",
                                         checks=versions.EQ(1, 4))

VERSION_15 = "1.5+"
_VERSION_15_CHECK = versions.Requirement(call=("AdapterRemoval", "--version"),
                                         search=r"ver. (\d+)\.(\d+)\.(\d+)",
                                         checks=versions.GE(1, 5, 0))


class SE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_files, output_prefix, output_format="bz2",
                  quality_offset=33, version=VERSION_15, dependencies=()):
        # See below for parameters in common between SE/PE
        cmd = _get_common_parameters(version)

        # Uncompressed reads (piped from 'paleomix cat')
        cmd.set_option("--file1",    "%(TEMP_IN_READS)s")

        # Prefix for output files
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        # Quality offset for Phred (or similar) scores
        assert quality_offset in (33, 64)
        cmd.set_option("--qualitybase", quality_offset)

        basename = os.path.basename(output_prefix)
        cmd.set_kwargs(# Only settings file is saved, rest is temporary files
                       OUT_SETTINGS        = output_prefix + ".settings",
                       TEMP_OUT_BASENAME   = basename,

                       # Named pipe for uncompressed input (paleomix cat)
                       TEMP_IN_READS       = "uncompressed_input",

                       # Named pipes for output of AdapterRemova
                       TEMP_OUT_LINK_1     = basename + ".truncated",
                       TEMP_OUT_LINK_2     = basename + ".discarded")

        return {"basename"      : basename,
                "version"       : version,
                "command"       : cmd}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._quality_offset = parameters.quality_offset
        self._basename = parameters.basename
        self._check_fastqs = _are_fastq_checks_required(parameters.version)

        zcat           = _build_cat_command(parameters.input_files, "uncompressed_input")
        zip_truncated  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".truncated")
        zip_discarded  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".discarded")
        adapterrm      = parameters.command.finalize()

        commands = ParallelCmds([adapterrm, zip_discarded, zip_truncated, zcat])
        CommandNode.__init__(self,
                             command      = commands,
                             description  = "<AdapterRM (SE): %s -> '%s.*'>" \
                                 % (fileutils.describe_files(parameters.input_files),
                                    parameters.output_prefix),
                             dependencies = parameters.dependencies)

    def _setup(self, config, temp):
        if self._check_fastqs:
            check_fastq_files(self.input_files, self._quality_offset)

        os.mkfifo(os.path.join(temp, self._basename + ".truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, "uncompressed_input"))

        CommandNode._setup(self, config, temp)


class PE_AdapterRemovalNode(CommandNode):
    @create_customizable_cli_parameters
    def customize(cls, input_files_1, input_files_2, output_prefix,
                  quality_offset=33, output_format="bz2", collapse=True,
                  version=VERSION_15, dependencies=()):
        cmd = _get_common_parameters(version)

        # Uncompressed mate 1 and 2 reads (piped from 'paleomix cat')
        cmd.set_option("--file1",    "%(TEMP_IN_READS_1)s")
        cmd.set_option("--file2",    "%(TEMP_IN_READS_2)s")

        # Prefix for output files, ensure that all end up in temp folder
        cmd.set_option("--basename", "%(TEMP_OUT_BASENAME)s")

        # Quality offset for Phred (or similar) scores
        assert quality_offset in (33, 64)
        cmd.set_option("--qualitybase", quality_offset)


        # Output files are explicity specified, to ensure that the order is
        # the same here as below. A difference in the order in which files are
        # opened can cause a deadlock, due to the use of named pipes
        # (see __init__).
        cmd.set_option("--output1", "%(TEMP_OUT_LINK_PAIR1)s")
        cmd.set_option("--output2", "%(TEMP_OUT_LINK_PAIR2)s")

        if collapse:
          cmd.set_option("--collapse")
          cmd.set_option("--outputcollapsed", "%(TEMP_OUT_LINK_ALN)s")
          if version == VERSION_15:
              cmd.set_option("--outputcollapsedtruncated",
                             "%(TEMP_OUT_LINK_ALN_TRUNC)s")

        cmd.set_option("--singleton", "%(TEMP_OUT_LINK_UNALN)s")
        cmd.set_option("--discarded", "%(TEMP_OUT_LINK_DISC)s")

        basename = os.path.basename(output_prefix)
        cmd.set_kwargs(# Only settings file is saved, rest is temporary files
                       OUT_SETTINGS=output_prefix + ".settings",
                       TEMP_OUT_BASENAME=basename,

                       # Named pipes for uncompressed input (paleomix cat)
                       TEMP_IN_READS_1="uncompressed_input_1",
                       TEMP_IN_READS_2="uncompressed_input_2",

                       # Named pipes for output of AdapterRemoval
                       TEMP_OUT_LINK_PAIR1=basename + ".pair1.truncated",
                       TEMP_OUT_LINK_PAIR2=basename + ".pair2.truncated",
                       TEMP_OUT_LINK_DISC=basename + ".discarded")

        if version == VERSION_15:
            cmd.set_kwargs(TEMP_OUT_LINK_ALN=basename + ".collapsed",
                           TEMP_OUT_LINK_ALN_TRUNC=basename + ".collapsed.truncated",
                           TEMP_OUT_LINK_UNALN=basename + ".singleton.truncated")
        elif version == VERSION_14:
            cmd.set_kwargs(TEMP_OUT_LINK_ALN=basename + ".singleton.aln.truncated",
                           TEMP_OUT_LINK_UNALN=basename + ".singleton.unaln.truncated")
        else:
            assert False

        return {"basename": basename,
                "command": cmd}

    @use_customizable_cli_parameters
    def __init__(self, parameters):
        self._quality_offset = parameters.quality_offset
        self._version = parameters.version
        self._basename = parameters.basename
        self._collapse = parameters.collapse
        self._check_fastqs = _are_fastq_checks_required(parameters.version)
        if len(parameters.input_files_1) != len(parameters.input_files_2):
            raise CmdError("Number of mate 1 files differ from mate 2 files: "
                           "%i != %i" % (len(parameters.input_files_1),
                                         len(parameters.input_files_2)))

        zcat_pair_1    = _build_cat_command(parameters.input_files_1, "uncompressed_input_1")
        zcat_pair_2    = _build_cat_command(parameters.input_files_2, "uncompressed_input_2")
        zip_pair_1     = _build_zip_command(parameters.output_format, parameters.output_prefix, ".pair1.truncated")
        zip_pair_2     = _build_zip_command(parameters.output_format, parameters.output_prefix, ".pair2.truncated")
        zip_discarded  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".discarded")
        adapterrm      = parameters.command.finalize()

        commands = [adapterrm, zip_pair_1, zip_pair_2]
        if parameters.version == VERSION_15:
            zip_unaligned  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".singleton.truncated")
            if parameters.collapse:
                zip_aln = _build_zip_command(parameters.output_format, parameters.output_prefix, ".collapsed")
                zip_aln_trunc = _build_zip_command(parameters.output_format, parameters.output_prefix, ".collapsed.truncated")
                commands += [zip_aln, zip_aln_trunc, zip_unaligned]
            else:
                commands += [zip_unaligned]
        else:
            zip_aln        = _build_zip_command(parameters.output_format, parameters.output_prefix, ".singleton.aln.truncated")
            zip_unaligned  = _build_zip_command(parameters.output_format, parameters.output_prefix, ".singleton.unaln.truncated")
            commands      += [zip_aln, zip_unaligned]
        commands += [zip_discarded, zcat_pair_1, zcat_pair_2]
        commands = ParallelCmds(commands)

        description  = "<AdapterRM (PE): %s -> '%s.*'>" \
            % (fileutils.describe_paired_files(parameters.input_files_1,
                                               parameters.input_files_2),
               parameters.output_prefix)

        CommandNode.__init__(self,
                             command=commands,
                             description=description,
                             dependencies=parameters.dependencies)

    def _setup(self, config, temp):
        if self._check_fastqs:
            check_fastq_files(self.input_files, self._quality_offset)

        os.mkfifo(os.path.join(temp, self._basename + ".discarded"))
        os.mkfifo(os.path.join(temp, self._basename + ".pair1.truncated"))
        os.mkfifo(os.path.join(temp, self._basename + ".pair2.truncated"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_1"))
        os.mkfifo(os.path.join(temp, "uncompressed_input_2"))

        if self._version == VERSION_15:
            if self._collapse:
               os.mkfifo(os.path.join(temp, self._basename + ".collapsed"))
               os.mkfifo(os.path.join(temp, self._basename + ".collapsed.truncated"))
            os.mkfifo(os.path.join(temp, self._basename + ".singleton.truncated"))
        else:
            os.mkfifo(os.path.join(temp, self._basename + ".singleton.aln.truncated"))
            os.mkfifo(os.path.join(temp, self._basename + ".singleton.unaln.truncated"))

        CommandNode._setup(self, config, temp)


def _build_cat_command(input_files, output_file):
    cat = factory.new("cat")
    cat.set_option("--output", "%(TEMP_OUT_CAT)s")
    cat.set_kwargs(TEMP_OUT_CAT=output_file)
    cat.add_multiple_values(input_files)

    return cat.finalize()


def _build_zip_command(output_format, prefix, name, output=None):
    if output_format not in ("gz", "bz2"):
        message = "Invalid output-format (%r), please select 'gz' or 'bz2'"
        raise CmdError(message % (output_format,))

    basename = os.path.basename(prefix)
    compress = factory.new("zip")
    compress.set_option("--format", output_format)
    compress.add_value("%(TEMP_IN_PIPE)s")
    compress.set_kwargs(TEMP_IN_PIPE=basename + name,
                        OUT_STDOUT=prefix + (output or name) + "." + output_format)

    return compress.finalize()


_DEPRECATION_WARNING_PRINTED = False


def _get_common_parameters(version):
    global _DEPRECATION_WARNING_PRINTED

    if version == VERSION_14:
        version_check = _VERSION_14_CHECK
    elif version == VERSION_15:
        version_check = _VERSION_15_CHECK
    else:
        raise CmdError("Unknown version: %s" % version)

    cmd = AtomicCmdBuilder("AdapterRemoval",
                           CHECK_VERSION=version_check)

    # Trim Ns at read ends
    cmd.set_option("--trimns", fixed=False)
    # Trim low quality scores
    cmd.set_option("--trimqualities", fixed=False)

    try:
        if not _DEPRECATION_WARNING_PRINTED and version_check.version < (2, 0):
            import pypeline.ui as ui
            ui.print_warn("\nWARNING: AdapterRemoval v1.5.x is deprecated;")
            ui.print_warn("         Upgrading to 2.1.x is strongly adviced!\n")
            ui.print_warn("         Download the newest version of AdapterRemoval at ")
            ui.print_warn("         https://github.com/MikkelSchubert/adapterremoval\n")

            _DEPRECATION_WARNING_PRINTED = True
    except versions.VersionRequirementError:
        pass

    return cmd


def _are_fastq_checks_required(version):
    if version == VERSION_14:
        return True
    elif version == VERSION_15:
        try:
            return _VERSION_15_CHECK.version < (2, 0)
        except versions.VersionRequirementError:
            return True
    else:
        raise CmdError("Unknown version: %s" % version)
