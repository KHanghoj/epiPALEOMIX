from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
import tempfile
import os


class GccorrectNode_good(CommandNode):
    def __init__(self, infile, temp_dir, dependencies=()):
        procs = list()
        rlmin, rlmax = infile.pop('MapMinMaxReadLength', [56, 57])
        bamfile = infile.pop("BamPath")
        output = tempfile.mkdtemp(dir=temp_dir, prefix='GCcorrect')
        destination_prefix = os.path.join(output,
                                          os.path.basename(bamfile)+'_')
        for rl in xrange(rlmin, rlmax+1, 1):
            destination = destination_prefix+str(rl)
            call = ['python', './gccorrect.py', '%(IN_BAM)s', '%(OUT_STDOUT)s']
            call.extend(("--ReadLength", rl))
            for item in infile.iteritems():
                call.extend(item)
            cmd = AtomicCmd(call, IN_BAM=bamfile, OUT_STDOUT=destination)
            procs.append(cmd)
        description = "<Gccorrect: '%s' -> '%s'>" % (bamfile,
                                                     (str(rlmin) + '_' +
                                                      str(rlmax)))
        CommandNode.__init__(self,
                             description=description,
                             command=ParallelCmds(procs),
                             dependencies=dependencies)


class GccorrectNode(CommandNode):
    def __init__(self, gc_dic, bamfile, destination_pref, rl, dependencies=()):
        destination = destination_pref+str(rl)
        call = ['python', './gccorrect.py', '%(IN_BAM)s', '%(OUT_STDOUT)s']
        call.extend(("--ReadLength", rl))
        for item in gc_dic.iteritems():
            call.extend(item)
        cmd = AtomicCmd(call, IN_BAM=bamfile, OUT_STDOUT=destination)
        description = "<Gccorrect: '%s' -> '%s'>" % (bamfile, str(rl))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class CreateGCModelNode(CommandNode):
    def __init__(self, source, pattern, destination_pref, dependencies=()):
        call = ['Rscript', './model_gc.R',
                '%(IN_SOURCE)s', str(pattern), '%(OUT_STDOUT)s']
        destination = os.path.join(destination_pref,
                                   'GC_Model_%s.txt' % (str(pattern)))
        cmd = AtomicCmd(call, IN_SOURCE=source, OUT_STDOUT=destination)
        description = "<CreateGCModel: '%s' -> '%s'" % (source,
                                                        destination)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)
