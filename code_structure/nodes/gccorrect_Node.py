from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
# from pypeline.atomiccmd.sets import ParallelCmds
import os


class GccorrectNode(CommandNode):
    def __init__(self, gc_dic, bamfile, destination_pref, rl, dependencies=()):
        destination = destination_pref+'_'+str(rl)
        call = ['python', './tools/gccorrect.py',
                '%(IN_BAM)s', '%(OUT_STDOUT)s']
        call.extend(("--ReadLength", rl))
        for option, argument in gc_dic.iteritems():
            if option.startswith('-'):
                call.extend((option, argument))
        cmd = AtomicCmd(call, IN_BAM=bamfile, OUT_STDOUT=destination)
        description = "<Gccorrect: '%s' -> '%s'>" % (bamfile, str(rl))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class CreateGCModelNode(CommandNode):
    def __init__(self, source_dest, outputfolder, dependencies=()):
        source, pattern = os.path.split(source_dest)
        call = ['Rscript', './tools/model_gc.R',
                '%(IN_SOURCE)s', str(pattern), outputfolder, '%(OUT_STDOUT)s']
        destination = os.path.join(outputfolder,
                                   'GC_Model_%s.txt' % (str(pattern)))
        cmd = AtomicCmd(call, IN_SOURCE=source, OUT_STDOUT=destination)
        description = "<CreateGCModel: '%s' -> '%s'" % (source,
                                                        destination)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)
