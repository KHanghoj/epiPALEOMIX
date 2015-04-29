from epipipe.node import CommandNode, Node
from epipipe.atomiccmd.command import AtomicCmd
from epipipe.atomiccmd.sets import ParallelCmds
import os


# class CleanFilesNodeold(CommandNode):
#     def __init__(self, mappa, unique, inbedfile, outbedfile, dependencies=()):
#         call = ['Rscript', './tools/cleanbedfiles.R',
#                 '%(IN_MAPPA)s', str(unique), '%(IN_BED)s', '%(OUT_BED)s']
#         cmd = AtomicCmd(call,
#                         IN_MAPPA=mappa,
#                         IN_BED=inbedfile,
#                         OUT_BED=outbedfile)
#         description = "<CreateGCModel: '%s' -> '%s', Uniqueness: '%s'" % \
#                       (inbedfile,
#                        outbedfile,
#                        unique)
#         CommandNode.__init__(self,
#                              description=description,
#                              command=cmd,
#                              dependencies=dependencies)


class CleanFilesNode(CommandNode):
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


import tools.splitbedfiles
                                  
class SplitBedFile(Node):
    def __init__(self, config, inbedfile, no_subbed, subnodes=(), dependencies=()):
        self.temp_root, self.infile = config.temp_root, inbedfile
        self.outputnames = self._createbednames(no_subbed)
        description = "<SplitBedFile: '%s' to '%s',"\
                      " Splitted in %s sub bed files" \
                      % (inbedfile, self.temp_root, str(no_subbed))
        Node.__init__(self,
                      description=description,
                      input_files=self.infile,
                      output_files=self.outputnames,
                      subnodes=subnodes,
                      dependencies=dependencies)
        assert isinstance(self.outputnames, list), "output has to be a list of str"

    def _run(self, _config, _temp):
        inputs = [self.infile]+self.outputnames
        tools.splitbedfiles.main(inputs)

    def _createbednames(self, no_subbed):
        ''' make each bedfile filename '''
        fmt, lst = '_0{}', []
        filena, fileext = os.path.splitext(os.path.basename(self.infile))
        for number in xrange(1,no_subbed+1):
            lst.append(os.path.join(self.temp_root, filena + fmt.format(number) + fileext))
        return lst


import tools.merge_datafiles
                                  
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
        tools.merge_datafiles.main(inputs)
