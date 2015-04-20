from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
# from pypeline.atomiccmd.sets import ParallelCmds
import os

GC_NAME = '_GCcorrect'


class GccorrectNode(CommandNode):
    def __init__(self, dat_bam, dat_make, rl, dependencies=()):
        dest = os.path.join(dat_bam.i_path,
                            dat_bam.bam_name+GC_NAME+'_'+str(rl))
        call = ['python', './tools/gccorrect.py', '%(IN_BAM)s',
                '%(OUT_STDOUT)s']
        call.extend(("--ReadLength", rl))
        for option, argument in dat_bam.retrievedat('GCcorrect').iteritems():
            if isinstance(option, str) and option.startswith('-'):
                call.extend((option, argument))
        cmd = AtomicCmd(call,
                        IN_BAM=dat_bam.baminfo["BamPath"],
                        OUT_STDOUT=dest)
        description = ("<Gccorrect: '%s' window length: '%s'>" %
                       (dat_bam.baminfo["BamPath"], str(rl)))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class CreateGCModelNode(CommandNode):
    def __init__(self, dat_bam, dependencies=(), subnodes=()):
        call = ['Rscript', './tools/model_gc.R',
                '%(IN_SOURCE)s', str(dat_bam.bam_name + GC_NAME),
                '%(OUT_FILEPATH)s', '%(OUT_PLOT)s']
        dest = os.path.join(dat_bam.o_path,
                            'GC_Model_%s.txt' % (str(dat_bam.bam_name)))
        plot_dest = os.path.splitext(dest)[0]+'.pdf'
        cmd = AtomicCmd(call,
                        IN_SOURCE=dat_bam.i_path,
                        OUT_FILEPATH=dest,
                        OUT_PLOT=plot_dest)
        description = "<CreateGCModel: '%s' -> '%s'" % (dat_bam.i_path,
                                                        dest)
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
