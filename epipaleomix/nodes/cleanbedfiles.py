#!/usr/bin/env python
from pypeline.common.fileutils import move_file, reroot_path
from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
import pypeline.common.versions as versions
from epipaleomix.tools import splitbedfiles
from epipaleomix.tools import merge_datafiles
import os, re

BEDTOOLS_VERSION = versions.Requirement(call   = ("bedtools", "--version"),
                                        search = r"bedtools v?(\d+)\.(\d+)\.(\d+)",
                                        checks = versions.GE(2, 15, 0))
PYTHON_VERSION = versions.Requirement(call   = ("python", "--version"),
                                        search = r"Python (\d+)\.(\d+)\.(\d+)",
                                        checks = versions.GE(2, 7, 0))


prefix = os.path.dirname(splitbedfiles.__file__)

class CleanFilesNode(CommandNode):
    def __init__(self, config, d_make, bedn, mappa, unique, dependencies=()):
        inbedfile = d_make.bedfiles[bedn]
        basename, extension = os.path.splitext(os.path.basename(inbedfile))
        bname = "{}_MappaOnly{}".format(basename, extension)
        outbedfile = os.path.join(config.temp_local, bname)
        d_make.bedfiles[bedn] = outbedfile
        
        ### outbedfile = os.path.join(config.temp_root, os.path.basename(inbedfile))

        call1 = ["python", os.path.join(prefix, 'filtermappa.py'),
                 "%(IN_MAPPA)s", str(unique)]
        call2 = ["bedtools", "intersect", "-wb", "-a", "stdin", "-b" ,"%(IN_BED)s"]
        ## call3 = ["sort",  "-V", "-k 4,4", "-k 5,5", "-k 2,2"]
        call3 = ["sort", "-k 4,4", "-k 5,5n", "-k 2,2n"]
        call4 = ["python", os.path.join(prefix, "updatebedcoord.py")]

        cmd1 = AtomicCmd(call1,
                         IN_MAPPA=mappa,
                         OUT_STDOUT=AtomicCmd.PIPE,
                         CHECK_VERSION = PYTHON_VERSION)
        cmd2 = AtomicCmd(call2,
                         IN_STDIN=cmd1,
                         IN_BED=inbedfile,
                         OUT_STDOUT=AtomicCmd.PIPE,
                         CHECK_VERSION = BEDTOOLS_VERSION)
        cmd3 = AtomicCmd(call3,
                         IN_STDIN=cmd2,
                         OUT_STDOUT=AtomicCmd.PIPE)
        cmd4 = AtomicCmd(call4,
                         IN_STDIN=cmd3,
                         OUT_STDOUT=outbedfile)

        paral_cmd = ParallelCmds([cmd1, cmd2, cmd3, cmd4])
        description = "<CLEANBEDFILES: '%s' -> '%s', Uniqueness: '%s'" % \
                      (inbedfile, outbedfile, unique)
        CommandNode.__init__(self,
                             description=description,
                             command=paral_cmd,
                             dependencies=dependencies)


class SplitBedFileNode(Node):
    def __init__(self, config, d_make, bedn, subnodes=(), dependencies=()):
        self.temp_local, self.inbedfile = config.temp_local, d_make.bedfiles[bedn]
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
                      subnodes=subnodes,
                      dependencies=dependencies)
        assert isinstance(self.outputnames, list), "output has to be a list of strings"

    def _run(self, _config, temp):
        temp_outputnames = [reroot_path(temp, dest) for dest in self.outputnames]
        inputs = [self.inbedfile]+temp_outputnames
        # inputs = [self.inbedfile]+self.outputnames
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
            lst.append(os.path.join(self.temp_local, filena + fmt.format(number) + fileext))
        return lst


class MergeDataFilesNode(Node):
    def __init__(self, d_bam, anal, bedn, subnodes=(), dependencies=()):
        self.infiles = [''.join(n.output_files) for n in subnodes]
        self.anal = anal
        analname = self._check_gccorr_name(anal, self.infiles[0])
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
                      subnodes=subnodes,
                      dependencies=dependencies)

    # def _run(self, _config, _temp):
    #     inputs = [self.anal, self.dest] + self.infiles  # INFILES IS A LIST ALREADY
    #     assert len(inputs) > 2, 'Need at least one output and one input'
    #     if self.anal == 'Phasogram':
    #         inputs.append('--merge')
    #     merge_datafiles.main(inputs)

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


    def _check_gccorr_name(self, anal, infile):
        curr_bamname, curr_anal, curr_bedname, _ = os.path.basename(infile).split('_')

        checkname = anal+"GCcorr"
        if re.search(checkname, curr_anal):
            anal += "GCcorr"
        assert anal == curr_anal, "Trying to merge different types of analyses: %s and %s" % (anal, curr_anal)
        return anal

