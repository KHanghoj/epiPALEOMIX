#!/usr/bin/env python
import os

import pypeline.common.versions as versions
from pypeline.common.fileutils import move_file, reroot_path
from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds

from epiomix.tools import splitbedfiles, \
    merge_datafiles


PYTHON_VERSION = versions.Requirement(call=("python", "--version"),
                                      search=r"Python (\d+)\.(\d+)\.(\d+)",
                                      checks=versions.GE(2, 7, 3))
PREFIX = os.path.dirname(splitbedfiles.__file__)


class CleanFilesNode(CommandNode):
    def __init__(self, config, d_make, bedn, mappa, unique, dependencies=()):
        inbedfile = d_make.bedfiles[bedn]
        basename, extension = os.path.splitext(os.path.basename(inbedfile))
        bname = "{}_MappaOnly{}".format(basename, extension)
        dest = os.path.join(config.temp_local, bname)

        d_make.bedfiles[bedn] = dest

        call1 = ["python", os.path.join(PREFIX, "intersectmappabed.py"),
                 "%(IN_BED)s", "%(IN_MAPPA)s", str(unique), "%(OUT_DEST)s"]

        
        cmd = AtomicCmd(call1,
                        IN_BED=inbedfile,
                        IN_MAPPA=mappa,
                        OUT_DEST=dest,
                        CHECK_VERSION=PYTHON_VERSION)

        description = ("<CLEANBEDFILES: '%s' -> '%s', Uniqueness: '%s'>" %
                       (inbedfile, dest, unique))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class SplitBedFileNode(Node):
    def __init__(self, config, d_make, bedn, dependencies=()):
        self.temp_local = config.temp_local
        self.inbedfile = d_make.bedfiles[bedn]
        self.no_subbed = config.max_threads
        self.outputnames = self._createbednames()
        d_make.bedfiles[bedn] = self.outputnames
        description = "<SplitBedFile: '%s' to '%s',"\
                      " Splitted in %s sub bed files" \
                      % (self.inbedfile, self.temp_local, str(self.no_subbed))
        Node.__init__(self,
                      description=description,
                      input_files=self.inbedfile,
                      output_files=self.outputnames,
                      dependencies=dependencies)
        assert isinstance(self.outputnames, list), \
            "output has to be a list of strings"

    def _run(self, _config, temp):
        temp_outputnames = [reroot_path(temp, dest)
                            for dest in self.outputnames]
        inputs = [self.inbedfile]+temp_outputnames
        splitbedfiles.main(inputs)

    def _teardown(self, _config, temp):
        for dest in self.outputnames:
            move_file(reroot_path(temp, dest), dest)

        Node._teardown(self, _config, temp)

    def _createbednames(self):
        ''' make each bedfile filename '''
        fmt, lst = '_0{}', []
        filena, fileext = os.path.splitext(os.path.basename(self.inbedfile))
        for number in xrange(0, self.no_subbed):
            lst.append(os.path.join(self.temp_local,
                                    filena + fmt.format(number) + fileext))
        return lst


class MergeDataFilesNode(Node):
    def __init__(self, d_bam, anal, bedn, dependencies=()):
        self.infiles = [''.join(n.output_files) for n in dependencies]
        self.anal = anal
        analname = self._check_gccorr_name(self.infiles[0])
        self.dest = os.path.join(d_bam.bam_output,
                                 d_bam.fmt.format(d_bam.bam_name,
                                                  analname,
                                                  bedn))
        assert self.infiles, "No temporary files to merge"
        if len(self.infiles) > 1:
            description = "<MergeDataFiles: '%s' ... '%s' ->  '%s'" % \
                          (os.path.basename(self.infiles[0]),
                           os.path.basename(self.infiles[-1]),
                           self.dest)
        else:
            description = "<MergeDataFiles: '%s' -> '%s'" % \
                          (os.path.basename(self.infiles[0]), self.dest)
        Node.__init__(self,
                      description=description,
                      input_files=self.infiles,
                      output_files=self.dest,
                      dependencies=dependencies)

    def _run(self, _config, temp):
        dest = reroot_path(temp, self.dest)
        inputs = [self.anal, dest] + self.infiles
        assert len(inputs) > 2, 'Need at least one output and one input'
        if self.anal == 'Phasogram':
            inputs.append('--merge')
        merge_datafiles.main(inputs)

    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self.dest), self.dest)
        Node._teardown(self, _config, temp)

    def _check_gccorr_name(self, infile):
        c_bamna, c_analname, c_bedna, _ = os.path.basename(infile).split('_')
        return c_analname
