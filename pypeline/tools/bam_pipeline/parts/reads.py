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

from pypeline.common.fileutils import missing_files
from pypeline.atomiccmd.builder import apply_options
from pypeline.nodes.adapterremoval import \
    SE_AdapterRemovalNode, \
    PE_AdapterRemovalNode, \
    VERSION_14, \
    VERSION_15
from pypeline.nodes.validation import \
    ValidateFASTQFilesNode


class Reads(object):
    def __init__(self, config, record, quality_offset):
        self.quality_offset = quality_offset
        self.files = {}
        self.stats = None
        self.nodes = ()

        tags = record["Tags"]
        self.folder = os.path.join(config.destination, tags["Target"], "reads",
                                   tags["SM"], tags["LB"], tags["PU_cur"])

        lane_type = record.get("Type")
        if lane_type == "Raw":
            self._init_raw_reads(record)
        elif lane_type == "Trimmed":
            self._init_pretrimmed_reads(record)
        else:
            assert False, "Unexpected data type in Reads(): %s" \
                % (repr(lane_type))

        for name in record["Options"]["ExcludeReads"]:
            self.files.pop(name, None)

        if config.allow_missing_input_files and self.nodes:
            input_missing = missing_files(self.nodes[0].input_files)
            output_missing = missing_files(self.nodes[0].output_files)
            if input_missing and not output_missing:
                self.nodes = ()

    def _init_pretrimmed_reads(self, record):
        self.files.update(record["Data"])
        output_file = os.path.join(self.folder, "reads.pretrimmed.validated")
        input_files = set()
        for (read_type, filename) in self.files.iteritems():
            if read_type == "Paired":
                input_files.add(filename.format(Pair=1))
                input_files.add(filename.format(Pair=2))
            else:
                input_files.add(filename)

        node = ValidateFASTQFilesNode(input_files=input_files,
                                      output_file=output_file,
                                      offset=self.quality_offset)
        self.nodes = (node,)

    def _init_raw_reads(self, record):
        # Support for older versions of the pipeline, which used ARv1.0 - 1.4
        version = VERSION_14
        if record["Options"]["AdapterRemoval"]["Version"] == "v1.5+":
            version = VERSION_15

        quality_offset = self.quality_offset
        if quality_offset == "Solexa":
            quality_offset = 64

        ar_options = dict(record["Options"]["AdapterRemoval"])
        # Setup of "--collapsed" is handled by the node itself
        collapse_reads = ar_options.pop("--collapse")
        collapse_reads = collapse_reads or collapse_reads is None

        init_args = {"output_prefix": os.path.join(self.folder, "reads"),
                     "output_format": record["Options"]["CompressionFormat"],
                     "quality_offset": quality_offset,
                     "version": version}
        output_tmpl = "{output_prefix}.%s.{output_format}".format(**init_args)

        if ("SE" in record["Data"]):
            self.files["Single"] = output_tmpl % ("truncated",)
            init_args["input_files"] = record["Data"]["SE"]
            command = SE_AdapterRemovalNode.customize(**init_args)
        else:
            if version is VERSION_14:
                self._set_adapterrm_v14_files(self.files, output_tmpl)
            else:
                self._set_adapterrm_v15_files(self.files, output_tmpl,
                                              collapse_reads)

            init_args["collapse"] = collapse_reads
            init_args["input_files_1"] = record["Data"]["PE_1"]
            init_args["input_files_2"] = record["Data"]["PE_2"]
            command = PE_AdapterRemovalNode.customize(**init_args)

        apply_options(command.command, ar_options)
        self.stats = os.path.join(self.folder, "reads.settings")
        self.nodes = (command.build_node(),)

    @classmethod
    def _set_adapterrm_v14_files(cls, files, output_tmpl):
        files["Single"] = output_tmpl % ("singleton.unaln.truncated",)
        files["Collapsed"] = output_tmpl % ("singleton.aln.truncated",)
        files["Paired"] = output_tmpl % ("pair{Pair}.truncated",)

    @classmethod
    def _set_adapterrm_v15_files(cls, files, output_tmpl, collapse_reads):
        files["Single"] = output_tmpl % ("singleton.truncated",)
        files["Paired"] = output_tmpl % ("pair{Pair}.truncated",)
        if collapse_reads:
            files["Collapsed"] = output_tmpl % ("collapsed",)
            files["CollapsedTruncated"] = output_tmpl % ("collapsed.truncated",)
