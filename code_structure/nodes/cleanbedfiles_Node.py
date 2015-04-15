from pypeline.node import CommandNode
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
