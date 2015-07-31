from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
import os
from epipaleomix.tools import splitbedfiles
from epipaleomix.tools import merge_datafiles

prefix = os.path.dirname(splitbedfiles.__file__)

class _CleanFilesNode_oldnode(CommandNode):
    def __init__(self, config, inbedfile, mappa, unique, dependencies=()):
        outbedfile = os.path.join(config.temp_root, os.path.basename(inbedfile))
        call1 = ["python", "./tools/filtermappa.py",
                 "%(IN_MAPPA)s", str(unique)]
        call2 = ["bedtools", "intersect", "-wa", "-u", "-a",
                 "%(IN_BED)s", "-b", "stdin"]
        cmd1 = AtomicCmd(call1,
                         IN_MAPPA=mappa,
                         OUT_STDOUT=AtomicCmd.PIPE)
        cmd2 = AtomicCmd(call2,
                         IN_STDIN=cmd1,
                         IN_BED=inbedfile,
                         OUT_STDOUT=outbedfile)
        paral_cmd = ParallelCmds([cmd1, cmd2])
        description = "<CLEANBEDFILES: '%s' -> '%s', Uniqueness: '%s'" % \
                      (inbedfile,
                       outbedfile,
                       unique)
        CommandNode.__init__(self,
                             description=description,
                             command=paral_cmd,
                             dependencies=dependencies)

class CleanFilesNode(CommandNode):
    def __init__(self, config, inbedfile, mappa, unique, dependencies=()):
        basename,extension = os.path.splitext(os.path.basename(inbedfile))
        bname = "{}_MappaOnly{}".format(basename, extension)
        outbedfile = os.path.join(config.temp_root, bname)
        ### outbedfile = os.path.join(config.temp_root, os.path.basename(inbedfile))
        call1 = ["python", os.path.join(prefix, 'filtermappa.py'),
                 "%(IN_MAPPA)s", str(unique)]
        call2 = ["bedtools", "intersect", "-wb", "-a", "stdin", "-b" ,"%(IN_BED)s"]
        call3 = ["sort",  "-V", "-k 4,4", "-k 5,5", "-k 2,2"]
        call4 = ["python", os.path.join(prefix, "updatebedcoord.py")]

        cmd1 = AtomicCmd(call1,
                         IN_MAPPA=mappa,
                         OUT_STDOUT=AtomicCmd.PIPE)
        cmd2 = AtomicCmd(call2,
                         IN_STDIN=cmd1,
                         IN_BED=inbedfile,
                         OUT_STDOUT=AtomicCmd.PIPE)
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

        


class SplitBedFile(Node):
    def __init__(self, config, inbedfile, bedn, subnodes=(), dependencies=()):
        self.temp_root, self.infile = config.temp_root, inbedfile[bedn]
        self.no_subbed = config.max_threads
        self.outputnames = self._createbednames()
        inbedfile[bedn] = self.outputnames
        description = "<SplitBedFile: '%s' to '%s',"\
                      " Splitted in %s sub bed files" \
                      % (self.infile, self.temp_root, str(self.no_subbed))
        Node.__init__(self,
                      description=description,
                      input_files=self.infile,
                      output_files=self.outputnames,
                      subnodes=subnodes,
                      dependencies=dependencies)
        assert isinstance(self.outputnames, list), "output has to be a list of strings"

    def _run(self, _config, _temp):
        inputs = [self.infile]+self.outputnames
        splitbedfiles.main(inputs)

    def _createbednames(self):
        ''' make each bedfile filename '''
        fmt, lst = '_0{}', []
        filena, fileext = os.path.splitext(os.path.basename(self.infile))
        for number in xrange(0, self.no_subbed):
            lst.append(os.path.join(self.temp_root, filena + fmt.format(number) + fileext))
        return lst


class MergeDataFiles(Node):
    def __init__(self, d_bam, anal, bedn, subnodes=(), dependencies=()):
        self.infiles = [''.join(n.output_files) for n in subnodes]
        self.anal = anal
        self.out = os.path.join(d_bam.o_path,
                                d_bam.fmt.format(d_bam.bam_name, anal, bedn))
        description = "<MergeDataFiles: '%s' to OUTPUT '%s'" % \
                      (self.infiles, self.out)
        Node.__init__(self,
                      description=description,
                      input_files=self.infiles,
                      output_files=self.out,
                      subnodes=subnodes,
                      dependencies=dependencies)

    def _run(self, _config, _temp):
        inputs = [self.anal, self.out] + self.infiles  # INFILES IS A LIST ALREADY
        assert len(inputs) > 2, 'Need at least one output and one input'
        if self.anal == 'Phasogram':
            inputs.append('--merge')
        merge_datafiles.main(inputs)
