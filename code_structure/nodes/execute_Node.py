from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
import os
#from copy import deepcopy
import tools.pileupdepth
import tools.nucleomap
import tools.methylmap
import tools.phasogram

        
MODULES = {
    'WriteDepth': tools.pileupdepth.main,
    'NucleoMap': tools.nucleomap.main,
    'MethylMap': tools.methylmap.main,
    'Phasogram': tools.phasogram.main
}



class GeneralExecuteNode_OLD(CommandNode):
    def __init__(self, aux_path, analysis_options, bedfile, dependencies=()):
        call = ['python', '%(AUX_PYTHON)s', '%(IN_BAM)s',
                '%(IN_BED)s', '%(OUT_FILEPATH)s']
        bam_name, anal_name, outfmt, options = deepcopy(analysis_options)
        bed_name, bed_path = bedfile
        base_bam = os.path.split(bam_name)[-1]
        destination = \
            os.path.join(bam_name, outfmt(base_bam, bed_name, anal_name))
        if not options.get('Apply_GC_Correction', True):
            options.pop('--GCmodel', None)

        for option, argument in options.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                call.extend((option, argument))
        cmd = AtomicCmd(call,
                        AUX_PYTHON=aux_path,
                        IN_BAM=options['BamPath'],
                        IN_BED=bed_path,
                        OUT_FILEPATH=destination+'.txt.gz')
        description = "<ANALYSIS:'%s', Bed:'%s', BAM: %s" % \
                      (anal_name, bed_name, base_bam)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class General_Plot_Node_OLD(CommandNode):
    def __init__(self, aux_path, infile, anal_name, dependencies=()):
        call = ['Rscript', '%(AUX_R)s', '%(IN_TXT)s', '%(OUT_FILEPATH)s']
        outfile = os.path.splitext(os.path.splitext(infile)[0])[0]+'.pdf'
        if anal_name in EXTRAS:
            call.extend(EXTRAS[anal_name])
        cmd = AtomicCmd(call,
                        AUX_R=aux_path,
                        IN_TXT=infile,
                        OUT_FILEPATH=outfile)
        description = "<PLOT_ANALYSIS:'%s', Infile:'%s', Outfile: %s" % \
                      (os.path.split(aux_path)[-1], infile, outfile)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class _GeneralExecuteNode(CommandNode):
    def __init__(self, anal, aux_python, d_bam, bed_name, bed_path,
                 dependencies=()):
        call = ['python', '%(AUX_PYTHON)s', '%(IN_BAM)s',
                '%(IN_BED)s', '%(OUT_FILEPATH)s']
        destination = \
            os.path.join(d_bam.o_path, d_bam.fmt.format(d_bam.bam_name,
                                                        anal, bed_name))

        opt_arg = d_bam.retrievedat(anal)
        gc_bool = opt_arg.get('Apply_GC_Correction', False)
        for option, argument in opt_arg.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                if option == '--GCmodel' and not gc_bool:
                    continue
                call.extend((option, argument))
        cmd = AtomicCmd(call,
                        AUX_PYTHON=aux_python,
                        IN_BAM=d_bam.baminfo['BamPath'],
                        IN_BED=bed_path,
                        OUT_FILEPATH=destination)
        description = "<ANALYSIS:'%s', BAM: %s, Bed:'%s'" % \
                      (anal, d_bam.bam_name, bed_name)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)



class GeneralExecuteNode1(CommandNode):
    '''this new node puts the file in temproray
    final file need to be a merge of the other files'''
    def __init__(self, anal, aux_python, d_bam, bed_name, bed_path,
                 dependencies=()):
        call = ['python', '%(AUX_PYTHON)s', '%(IN_BAM)s',
                '%(IN_BED)s', '%(OUT_FILEPATH)s']
        destination = \
            os.path.join(d_bam.i_path, d_bam.fmt.format(d_bam.bam_name,
                                                        anal, bed_name))
        # d_bam.o_path
        opt_arg = d_bam.retrievedat(anal)
        gc_bool = opt_arg.get('Apply_GC_Correction', False)
        for option, argument in opt_arg.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                if option == '--GCmodel' and not gc_bool:
                    continue
                call.extend((option, argument))
        cmd = AtomicCmd(call,
                        AUX_PYTHON=aux_python,
                        IN_BAM=d_bam.baminfo['BamPath'],
                        IN_BED=bed_path,
                        OUT_FILEPATH=destination)
        description = "<ANALYSIS:'%s', BAM: %s, Bed:'%s'" % \
                      (anal, d_bam.bam_name, bed_name)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class GeneralExecuteNode(Node):
    '''this new node puts the file in temporary
    final file need to be a merge of the other files'''
    def __init__(self, anal, d_bam, bed_name, bed_path, dependencies=()):
        self.analysis = MODULES[anal]
        self.infile, self.d_bam = d_bam.baminfo['BamPath'], d_bam
        self.dest = os.path.join(d_bam.i_path,
                                 d_bam.fmt.format(d_bam.bam_name, anal, bed_name))
        self.inputs = [self.infile, bed_path, self.dest]
        self._add_options(anal)
        description = "<ANALYSIS:'%s', BAM: %s, Bed:'%s'" % \
                      (anal, d_bam.bam_name, bed_name)
        Node.__init__(self,
                      description=description,
                      input_files=[self.infile, bed_path],
                      output_files=self.dest,
                      dependencies=dependencies)

    def _run(self, _config, _temp):
        self.analysis(self.inputs)

    def _add_options(self, name):
        opt_arg = self.d_bam.retrievedat(name)
        gc_bool = opt_arg.get('Apply_GC_Correction', False)
        for option, argument in opt_arg.items():
            if isinstance(option, str) and option.startswith('-'):
                if option == '--GCmodel' and not gc_bool:
                    continue
                if not isinstance(argument, str):
                    argument = str(argument)
                self.inputs.extend((option, argument))


EXTRAS = {'Phasogram': ['1000']}

RPATHS = {'Phasogram': './tools/phaso.R',
          'WriteDepth': '_',
          'NucleoMap': './tools/nucleo_merge.R',
          'MethylMap': './tools/methyl_merge.R'
}



class General_Plot_Node(CommandNode):
    def __init__(self, infile, anal_name, dependencies=()):
        call = ['Rscript', '%(AUX_R)s', '%(IN_TXT)s', '%(OUT_FILEPATH)s']
        outfile = os.path.splitext(os.path.splitext(infile)[0])[0]+'.pdf'
        r_anal = RPATHS[anal_name]
        if anal_name in EXTRAS:
            call.extend(EXTRAS[anal_name])
        cmd = AtomicCmd(call,
                        AUX_R=r_anal,
                        IN_TXT=infile,
                        OUT_FILEPATH=outfile)
        description = "<PLOT_ANALYSIS:'%s', Infile:'%s', Outfile: %s" % \
                      (os.path.split(r_anal)[-1], infile, outfile)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)
