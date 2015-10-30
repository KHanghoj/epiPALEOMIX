import os
from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
from epipaleomix.tools import splitbedfiles
from epipaleomix.tools import merge_datafiles

import pypeline.common.versions as versions


BEDTOOLS_VERSION = versions.Requirement(call   = ("bedtools", "--version"),
                                        search = r"bedtools v?(\d+)\.(\d+)\.(\d+)",
                                        checks = versions.GE(2, 15, 0))
PYTHON_VERSION = versions.Requirement(call   = ("python", "--version"),
                                        search = r"Python (\d+)\.(\d+)\.(\d+)",
                                        checks = versions.GE(2, 7))


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
        call3 = ["sort",  "-V", "-k 4,4", "-k 5,5", "-k 2,2"]
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

    def _run(self, _config, _temp):
        inputs = [self.inbedfile]+self.outputnames
        splitbedfiles.main(inputs)

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
        opt_arg = d_bam.retrievedat(anal)
        self.analname = anal + 'GCcorr' if opt_arg.get('Apply_GC_Correction', False) else anal
        self.dest = os.path.join(d_bam.bam_output,
                                d_bam.fmt.format(d_bam.bam_name,
                                                 self.analname,
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

    def _run(self, _config, _temp):
        inputs = [self.anal, self.dest] + self.infiles  # INFILES IS A LIST ALREADY
        assert len(inputs) > 2, 'Need at least one output and one input'
        if self.anal == 'Phasogram':
            inputs.append('--merge')
        merge_datafiles.main(inputs)
