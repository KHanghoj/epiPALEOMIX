from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds

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
    def __init__(self, mappa, unique, inbedfile, outbedfile, dependencies=()):
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
    def __init__(self, temp_root, inbedfile, no_subbed=3, subnodes=(), dependencies=()):
        self.temp_root, self.infile = temp_root, inbedfile
        self.outputnames = self._createbednames(no_subbed)
        description = "<SplitBedile: '%s' to temproot '%s'" % \
                      (inbedfile,
                       temp_root)
        Node.__init__(self,
                      description=description,
                      input_files=self.infile,
                      output_files=self.outputnames,
                      subnodes=subnodes
                      dependencies=dependencies)

    def _run(self, _config, _temp):
        tools.splitbedfiles.main(self.infile, self.outputnames)

    def _createbednames(self, no_subbed):
        ''' make each bedfile filename '''
        fmt, lst = '_0{}', []
        filena, fileext = os.path.splitext(os.path.basename(self.infile))
        for number in xrange(1,no_subbed+1):
            lst.append(os.path.join(self.temp_root, filena + fmt.format(number) + fileext))
        return lst
