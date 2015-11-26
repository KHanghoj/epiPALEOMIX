#!/usr/bin/env python
from pypeline.node import Node
from pypeline.common.fileutils import move_file, reroot_path
import os

from epiomix.tools import \
    pileupdepth, \
    nucleomap, \
    methylmap, \
    phasogram

# EXTRAS = {'Phasogram': ['1000']}

# RPATHS = {'Phasogram': 'phaso.R',
#           'WriteDepth': '_',
#           'NucleoMap': 'nucleo_merge.R',
#           'MethylMap': 'methyl_merge.R'
# }

PREFIX = os.path.dirname(nucleomap.__file__)
        
MODULES = {
    'WriteDepth': pileupdepth.main,
    'NucleoMap': nucleomap.main,
    'MethylMap': methylmap.main,
    'Phasogram': phasogram.main
}


class GeneralExecuteNode(Node):
    '''this new node puts the file in temporary
    final file need to be a merge of the other files'''
    def __init__(self, anal, d_bam, bed_name, bed_path, gcnode, splitnode,  dependencies=()):
        self.analysis = MODULES[anal]
        self.infile, self.d_bam = d_bam.baminfo['BamPath'], d_bam
        self.bed_path, self.anal = bed_path, anal
        dependencies, analname = self._correct_subnodes(gcnode, splitnode, anal)

        out_f_name = d_bam.fmt.format(d_bam.bam_name, analname, bed_name)
        self.dest = os.path.join(d_bam.bam_temp_local, out_f_name)

        description = "<ANALYSIS:'%s', BAM: %s, Bed:'%s'" % \
                      (analname, d_bam.bam_name, bed_name)
        Node.__init__(self,
                      description=description,
                      input_files=[self.infile, bed_path],
                      output_files=self.dest,
                      dependencies=dependencies)

    def _run(self, _config, temp):
        dest = reroot_path(temp, self.dest)
        self.inputs = [self.infile, self.bed_path, dest]
        self._add_options(self.anal)
        self.analysis(self.inputs)

    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self.dest), self.dest)
        Node._teardown(self, _config, temp)

    def _add_options(self, name):
        opt_arg = self.d_bam.retrievedat(name)
        gc_bool = opt_arg.get('Apply_GC_Correction', False)
        for option, argument in opt_arg.items():
            if isinstance(option, str) and option.startswith('-'):
                if option == '--GCmodel' and not gc_bool:
                    continue
                if not isinstance(argument, str):
                    argument = str(argument)
                self.inputs.extend([option, argument])

    def _correct_subnodes(self, gcnode, splitnode, name):
        ''' As not all analyses requires to wait for gccorrection model.
            We remove it from the tuple of dependencies'''
        opt_arg = self.d_bam.retrievedat(name)
        if opt_arg.get('Apply_GC_Correction', False) and gcnode:
            return gcnode+splitnode, name+'GCcorr'
        else:
            return splitnode, name


# class GeneralPlotNode(CommandNode):
#     def __init__(self, infile, anal_name, dependencies=()):
#         call = ['Rscript', '%(AUX_R)s', str(50), '%(IN_TXT)s', '%(OUT_FILEPATH)s']
#         outfile = os.path.splitext(os.path.splitext(infile)[0])[0]+'.pdf'
#         r_anal = os.path.join(PREFIX, RPATHS[anal_name])
#         if anal_name in EXTRAS:
#             call.extend(EXTRAS[anal_name])
#         cmd = AtomicCmd(call,
#                         AUX_R=r_anal,
#                         IN_TXT=infile,
#                         OUT_FILEPATH=outfile)
#         description = "<PLOT_ANALYSIS:'%s', Infile:'%s', Outfile: %s" % \
#                       (os.path.split(r_anal)[-1], infile, outfile)
#         CommandNode.__init__(self,
#                              description=description,
#                              command=cmd,
#                              dependencies=dependencies)
