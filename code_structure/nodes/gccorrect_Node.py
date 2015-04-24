from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
# from pypeline.atomiccmd.builder import AtomicCmdBuilder
import os
import tools.gccorrect
GC_NAME = '_GCcorrect'


class GccorrectNode_OLD(CommandNode):
    def __init__(self, d_bam, rl, offset=None, subnodes=(), dependencies=()):
        dest = os.path.join(d_bam.i_path,
                            d_bam.bam_name+GC_NAME+'_'+str(rl))
        call = ['python', './tools/gccorrect.py', '%(IN_BAM)s',
                '%(OUT_STDOUT)s']
        if offset:
            call.extend(("--OffSet", list(offset.output_files)[0]))
            dest += '_finescale'
        call.extend(("--ReadLength", rl))
        for option, argument in d_bam.retrievedat('GCcorrect').iteritems():
            if isinstance(option, str) and option.startswith('-'):
                call.extend((option, argument))
        cmd = AtomicCmd(call,
                        IN_BAM=d_bam.baminfo["BamPath"],
                        OUT_STDOUT=dest)
        description = ("<Gccorrect: '%s' window length: '%s'>" %
                       (d_bam.baminfo["BamPath"], str(rl)))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             subnodes=subnodes,
                             dependencies=dependencies)


class _GccorrectNode(CommandNode):
    def __init__(self, d_bam, rl, subnodes=(), dependencies=()):
        dest = os.path.join(d_bam.i_path,
                            d_bam.bam_name+GC_NAME+'_'+str(rl))
        call = ['python', './tools/gccorrect.py', '%(IN_BAM)s',
                '%(OUT_STDOUT)s']
        if subnodes:
            call.extend(("--OffSet", list(subnodes[0].output_files)[0]))
            dest += '_finescale'
        call.extend(("--ReadLength", rl))
        for option, argument in d_bam.retrievedat('GCcorrect').iteritems():
            if isinstance(option, str) and option.startswith('-'):
                call.extend((option, argument))
        cmd = AtomicCmd(call,
                        IN_BAM=d_bam.baminfo["BamPath"],
                        OUT_STDOUT=dest)
        description = ("<Gccorrect: '%s' window length: '%s'>" %
                       (d_bam.baminfo["BamPath"], str(rl)))
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             subnodes=subnodes,
                             dependencies=dependencies)


class GccorrectNode(Node):
    def __init__(self, d_bam, rl, subnodes=(), dependencies=()):
        self.defrozen_out = os.path.join(d_bam.i_path,
                                         d_bam.bam_name+GC_NAME+'_'+str(rl))
        self.rl, self.d_bam, self.subns = rl, d_bam, subnodes
        if self.subns:
            self.defrozen_out += '_finescale'
        description = ("<Gccorrect: '%s' window length: '%s'>" %
                       (self.defrozen_out, rl))

        Node.__init__(self,
                      description=description,
                      input_files=self.d_bam.baminfo["BamPath"],
                      output_files=self.defrozen_out,
                      subnodes=subnodes,
                      dependencies=dependencies)
        assert len(self.output_files) == 1, self.output_files

    def _run(self, _config, _temp):
        self.inputs = [self.d_bam.baminfo["BamPath"], self.defrozen_out]
        if self.subns:
            offsetfile = (''.join(node.output_files) for node in self.subns)
            self.inputs.extend(("--OffSet", str(offsetfile.next())))
        self._add_options('GCcorrect')
        self.inputs.extend(("--ReadLength", str(self.rl)))
        tools.gccorrect.main(self.inputs)

    def _add_options(self, name):
        optargs = self.d_bam.retrievedat(name)
        for option, argument in optargs.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                if not isinstance(argument, str):
                    argument = str(argument)
                self.inputs.extend((option, argument))


# class GccorrectNode_Mid(CommandNode):
#     def __init__(self, d_bam, subnodes=()):
#         infiles = [list(node.output_files)[0] for node in subnodes]
#         dest = os.path.join(d_bam.i_path,
#                             d_bam.bam_name+'_MID'+GC_NAME+'.txt')
#         call = ['python', './tools/gccorrect_mid.py',
#                 '%(OUT_FILEPATH)s']
#         cmd = AtomicCmdBuilder(call, OUT_FILEPATH=dest)
#         cmd.add_multiple_values(infiles)
#         description = ("<GccorrectMid: 'infiles' to: '%s'>" %
#                        (dest))
#         CommandNode.__init__(self,
#                              description=description,
#                              command=cmd.finalize(),
#                              subnodes=subnodes,
#                              dependencies=())


class CreateGCModelNode(CommandNode):
    def __init__(self, d_bam, subnodes=()):
        call = ['Rscript', './tools/model_gc.R',
                '%(IN_SOURCE)s', str(d_bam.bam_name+GC_NAME),
                '%(OUT_FILEPATH)s', '%(OUT_PLOT)s']
        dest = os.path.join(d_bam.o_path,
                            'GC_Model_%s.txt' % (str(d_bam.bam_name)))
        plot_dest = os.path.splitext(dest)[0]+'.pdf'
        cmd = AtomicCmd(call,
                        IN_SOURCE=d_bam.i_path,
                        OUT_FILEPATH=dest,
                        OUT_PLOT=plot_dest)
        description = "<CreateGCModel: '%s' -> '%s'" % (d_bam.i_path,
                                                        dest)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=(),
                             subnodes=subnodes)


import tools.gccorrect_mid


class GccorrectNode_Mid(Node):
    def __init__(self, d_bam, subnodes=()):
        dest = os.path.join(d_bam.i_path,
                            d_bam.bam_name+'_MID'+GC_NAME+'.txt')
        description = ("<GccorrectMid: 'infiles' to: '%s'>" %
                       (dest))
        self.defrozen_in = [''.join(node.output_files) for node in subnodes]
        # self.defrozen_out = dest
        self.defrozen_out = [dest]
        Node.__init__(self,
                      description=description,
                      input_files=self.defrozen_in,
                      output_files=self.defrozen_out,
                      subnodes=subnodes,
                      dependencies=())
        assert len(self.output_files) == 1, self.output_files

    def _run(self, _config, _temp):
        tools.gccorrect_mid.main(self.defrozen_out+self.defrozen_in)

# class MyNode(Node):
#     def __init__(self, pattern, io_paths, dependencies=(), subnodes=()):
#         pass

#     def _run(self):
#         # this method overrides the empty Node _run.
#         # has to return the out of the imported script/python module.
