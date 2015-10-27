from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
# from pypeline.atomiccmd.builder import AtomicCmdBuilder
import os
from epipaleomix.tools import gccorrect_mid, gccorrect
GC_NAME = '_GCcorrect'


class GccorrectNode(Node):
    def __init__(self, d_bam, rl, subnodes=(), dependencies=()):
        self.defrozen_out = os.path.join(d_bam.i_path,
                                         d_bam.bam_name+GC_NAME+'_'+str(rl))
        self.rl, self.d_bam, self.subns = rl, d_bam, subnodes
        firstchrom, secondchrom = self.d_bam.opts['GCcorrect']['ChromUsed']
        if self.subns:  ## finetune step run if subnodes
            self.defrozen_out += '_finescale'
            self.chromtobeanalyzed = secondchrom
        else:
            self.chromtobeanalyzed = firstchrom
        description = ("<Gccorrect: '%s' window length: '%s' based on chromosome %s>" %
                       (self.defrozen_out, rl, self.chromtobeanalyzed))

        Node.__init__(self,
                      description=description,
                      input_files=self.d_bam.baminfo["BamPath"],
                      output_files=self.defrozen_out,
                      subnodes=subnodes,
                      dependencies=dependencies)
        assert len(self.output_files) == 1, self.output_files

    def _run(self, _config, _temp):

        self.inputs = [self.d_bam.baminfo["BamPath"], self.defrozen_out]
        if self.subns:  ## finetune step run if subnodes
            offsetfile = (''.join(node.output_files) for node in self.subns)
            self.inputs.extend(("--OffSet", str(offsetfile.next())))

        self.inputs.extend(("--ChromUsed", self.chromtobeanalyzed))
        self._add_options('GCcorrect')
        self.inputs.extend(("--ReadLength", str(self.rl)))
        gccorrect.main(self.inputs)

    def _add_options(self, name):
        optargs = self.d_bam.retrievedat(name)
        for option, argument in optargs.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                if not isinstance(argument, str):
                    argument = str(argument)
                self.inputs.extend((option, argument))
        

class CreateGCModelNode(CommandNode):
    def __init__(self, d_bam, subnodes=()):
        aux_r = os.path.join(os.path.dirname(gccorrect.__file__),
                             'model_gc.R')
        call = ['Rscript', aux_r, '%(IN_SOURCE)s', str(d_bam.bam_name+GC_NAME),
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
        d_bam.opts['BamInfo']['--GCmodel'] = dest


class GccorrectMidNode(Node):
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
        gccorrect_mid.main(self.defrozen_out+self.defrozen_in)
