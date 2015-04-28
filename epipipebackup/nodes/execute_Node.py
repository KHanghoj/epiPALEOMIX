from epipipe.node import CommandNode, Node
from epipipe.atomiccmd.command import AtomicCmd
import os
from epipipe.tools import \
    pileupdepth, \
    nucleomap, \
    methylmap, \
    phasogram

EXTRAS = {'Phasogram': ['1000']}

RPATHS = {'Phasogram': 'phaso.R',
          'WriteDepth': '_',
          'NucleoMap': 'nucleo_merge.R',
          'MethylMap': 'methyl_merge.R'
}

prefix = os.path.dirname(Nucleomap.__file__)
        
MODULES = {
    'WriteDepth': pileupdepth.main,
    'NucleoMap': nucleomap.main,
    'MethylMap': methylmap.main,
    'Phasogram': phasogram.main
}


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


class General_Plot_Node(CommandNode):
    def __init__(self, infile, anal_name, dependencies=()):
        call = ['Rscript', '%(AUX_R)s', '%(IN_TXT)s', '%(OUT_FILEPATH)s']
        outfile = os.path.splitext(os.path.splitext(infile)[0])[0]+'.pdf'
        r_anal = os.path.join(prefix, RPATHS[anal_name])
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
