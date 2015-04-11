from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd


class CleanFilesNode(CommandNode):
    def __init__(self, mappa, unique, inbedfile, outbedfile, dependencies=()):
        call = ['Rscript', './tools/cleanbedfiles.R',
                '%(IN_MAPPA)s', str(unique), '%(IN_BED)s', '%(OUT_BED)s']
        cmd = AtomicCmd(call, IN_MAPPA=mappa,
                        IN_BED=inbedfile,  OUT_BED=outbedfile)
        description = "<CreateGCModel: '%s' -> '%s', Uniqueness: '%s'" % \
                      (inbedfile,
                       outbedfile,
                       unique)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)
