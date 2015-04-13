from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
# from pypeline.atomiccmd.sets import ParallelCmds
import os


class GccorrectNode(CommandNode):
    def __init__(self, gc_dic, destination_pref, rl, dependencies=()):
        destination = destination_pref+'_'+str(rl)
        call = ['python', './tools/gccorrect.py', '%(IN_BAM)s',
                '%(OUT_STDOUT)s']
        call.extend(("--ReadLength", rl))
        for option, argument in gc_dic.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                call.extend((option, argument))
        cmd = AtomicCmd(call, IN_BAM=gc_dic["BamPath"],
                        OUT_STDOUT=destination)
        description = "<Gccorrect: '%s' -> '%s'>" % (gc_dic["BamPath"], str(rl))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class CreateGCModelNode(CommandNode):
    def __init__(self, pattern, io_paths, dependencies=()):
        call = ['Rscript', './tools/model_gc.R',
                '%(IN_SOURCE)s', str(pattern),
                io_paths['o_out'], '%(OUT_STDOUT)s']
        destination = os.path.join(io_paths['o_out'],
                                   'GC_Model_%s.txt' % (str(pattern)))
        cmd = AtomicCmd(call, IN_SOURCE=io_paths['temp'],
                        OUT_STDOUT=destination)
        description = "<CreateGCModel: '%s' -> '%s'" % (io_paths['temp'],
                                                        destination)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)
