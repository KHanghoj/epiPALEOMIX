from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
# from pypeline.atomiccmd.sets import ParallelCmds
import os


class GccorrectNode(CommandNode):
    def __init__(self, gc_dic, suffix, io_paths, rl, dependencies=()):
        destination = os.path.join(io_paths['temp'], suffix+'_'+str(rl))
        call = ['python', './tools/gccorrect.py', '%(IN_BAM)s',
                '%(OUT_STDOUT)s']
        call.extend(("--ReadLength", rl))
        for option, argument in gc_dic.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                call.extend((option, argument))
        cmd = AtomicCmd(call,
                        IN_BAM=gc_dic["BamPath"],
                        OUT_STDOUT=destination)
        description = "<Gccorrect: '%s' -> '%s'>" % (gc_dic["BamPath"], str(rl))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class CreateGCModelNode(CommandNode):
    def __init__(self, pattern, io_paths, dependencies=(), subnodes=()):
        call = ['Rscript', './tools/model_gc.R',
                '%(IN_SOURCE)s', str(pattern),
                '%(OUT_FILEPATH)s', '%(OUT_PLOT)s']
        destination = os.path.join(io_paths['o_out'],
                                   'GC_Model_%s.txt' % (str(pattern)))
        plot_dest = os.path.splitext(destination)[0]+'.pdf'
        cmd = AtomicCmd(call,
                        IN_SOURCE=io_paths['temp'],
                        OUT_FILEPATH=destination,
                        OUT_PLOT=plot_dest)
        description = "<CreateGCModel: '%s' -> '%s'" % (io_paths['temp'],
                                                        destination)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies,
                             subnodes=subnodes)


class CreateGCModelNode_old(CommandNode):
    def __init__(self, pattern, io_paths, dependencies=(), subnodes=()):
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
                             dependencies=dependencies,
                             subnodes=subnodes)


# class MyNode(Node):
#     def __init__(self, pattern, io_paths, dependencies=(), subnodes=()):
#         pass

#     def _run(self):
#         # this method overrides the empy Node _run.
#         # has to return the out of the imported script/python module.
