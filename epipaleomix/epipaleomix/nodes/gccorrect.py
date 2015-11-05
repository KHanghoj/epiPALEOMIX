import os
from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
import pypeline.common.versions as versions

from epipaleomix.tools import gccorrect_mid, gccorrect

Rscript_VERSION = versions.Requirement(call   = ("Rscript", "--version"),
                                       search = r"version (\d+)\.(\d+)\.(\d+)",
                                       checks = versions.GE(2, 15, 3))

GC_NAME = '_GCcorrect'


class GccorrectNode(Node):
    def __init__(self, d_bam, rl, subnodes=(), dependencies=()):
        self.dest = os.path.join(d_bam.bam_temp_local,
                                 d_bam.bam_name+GC_NAME+'_'+str(rl))
        self.rl, self.d_bam, self.subns = rl, d_bam, subnodes
        if self.subns:  ## finetune step run if subnodes
            self.dest += '_finescale'
        description = ("<Gccorrect: '%s' window length: '%s' based on chromosome %s>" %
                       (self.dest, rl, self.d_bam.opts['GCcorrect']['--ChromUsed']))

        Node.__init__(self,
                      description=description,
                      input_files=self.d_bam.baminfo["BamPath"],
                      output_files=self.dest,
                      subnodes=subnodes,
                      dependencies=dependencies)
        assert len(self.output_files) == 1, self.output_files


    def _run(self, _config, _temp):
        self.inputs = [self.d_bam.baminfo["BamPath"], self.dest]
        if self.subns:  ## finetune step run if subnodes
            offsetfile = (''.join(node.output_files) for node in self.subns)
            self.inputs.extend(("--OffSet", str(offsetfile.next())))

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
        # aux_r = os.path.join(os.path.dirname(gccorrect.__file__),
        #                      'model_gc_individualreadlength.R')
        call = ['Rscript', aux_r, '%(IN_SOURCE)s', str(d_bam.bam_name+GC_NAME),
                '%(OUT_FILEPATH)s', '%(OUT_PLOT)s']
        dest = os.path.join(d_bam.bam_output,
                            '%s_GC_Model.txt' % (str(d_bam.bam_name)))
        plot_dest = os.path.splitext(dest)[0]+'.pdf'
        cmd = AtomicCmd(call,
                        IN_SOURCE=d_bam.bam_temp_local,
                        OUT_FILEPATH=dest,
                        OUT_PLOT=plot_dest,
                        CHECK_VERSION = Rscript_VERSION)
        d_bam.opts['BamInfo']['--GCmodel'] = dest
        description = "<CreateGCModel: '%s' -> '%s'" % (d_bam.bam_temp_local,
                                                        dest)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=(),
                             subnodes=subnodes)


class GccorrectMidNode(Node):
    def __init__(self, d_bam, subnodes=()):
        dest = os.path.join(d_bam.bam_temp_local,
                            d_bam.bam_name+'_MID'+GC_NAME+'.txt')
        description = ("<GccorrectMid: 'infiles' to: '%s'>" %
                       (dest))
        self.infiles = [''.join(node.output_files) for node in subnodes]
        self.dest = [dest]
        assert len(self.dest) == 1, self.dest 
        Node.__init__(self,
                      description=description,
                      input_files=self.infiles,
                      output_files=self.dest,
                      subnodes=subnodes,
                      dependencies=())

    def _run(self, _config, _temp):
        gccorrect_mid.main(self.dest+self.infiles)


### this is for individual readlength gccorrection
class _GccorrectNode(Node):
    def __init__(self, d_bam, rl, subnodes=(), dependencies=()):
        self.dest = os.path.join(d_bam.bam_temp_local,
                                 d_bam.bam_name+GC_NAME+'_'+str(rl))
        self.rl, self.d_bam, self.subns = rl, d_bam, subnodes
        description = ("<Gccorrect: '%s' window length: '%s' based on chromosome %s>" %
                       (self.dest, rl, self.d_bam.opts['GCcorrect']['--ChromUsed']))

        Node.__init__(self,
                      description=description,
                      input_files=self.d_bam.baminfo["BamPath"],
                      output_files=self.dest,
                      subnodes=subnodes,
                      dependencies=dependencies)
        assert len(self.output_files) == 1, self.output_files


    def _run(self, _config, _temp):
        self.inputs = [self.d_bam.baminfo["BamPath"], self.dest]
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


