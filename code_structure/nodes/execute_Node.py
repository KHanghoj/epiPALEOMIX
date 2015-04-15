from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
import os
from copy import deepcopy


class GeneralExecuteNode(CommandNode):
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


EXTRAS = {'Phasogram': ['1000']}
# 'WriteDepth': ['./tools/pileupdepth.py', '_'],
# 'NucleoMap': ['./tools/nucleomap.py', './tools/nucleo_merge.R'],
# 'MethylMap': ['./tools/methylmap.py', './tools/methyl_merge.R']
# }


class General_Plot_Node(CommandNode):
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
