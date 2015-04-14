from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
# from pypeline.atomiccmd.sets import ParallelCmds
import os
from copy import deepcopy


class PhasoNode(CommandNode):
    def __init__(self, call, analysis_options, bedfile, dependencies=()):
        self.call = call
        bam_name, anal_name, outfmt, options = deepcopy(analysis_options)
        bed_name, bed_path = bedfile
        base_bam = os.path.split(bam_name)[-1]
        destination = \
            os.path.join(bam_name, outfmt(base_bam, bed_name, anal_name))
        if not options['Apply_GC_Correction']:
            options.pop('--GCmodel', None)

        for option, argument in options.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                self.call.extend((option, argument))
        cmd = AtomicCmd(self.call, IN_BAM=options['BamPath'],
                        IN_BED=bed_path,
                        OUT_FILEPATH=destination+'.txt.gz')
        description = "<ANALYSIS:'%s', Bed:'%s', BAM: %s" % \
                      (anal_name, bed_name, base_bam)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class WriteNode(CommandNode):
    def __init__(self, call, analysis_options, bedfile, dependencies=()):
        self.call = call
        bam_name, anal_name, outfmt, options = deepcopy(analysis_options)
        bed_name, bed_path = bedfile
        base_bam = os.path.split(bam_name)[-1]
        destination = \
            os.path.join(bam_name, outfmt(base_bam, bed_name, anal_name))
        if not options['Apply_GC_Correction']:
            options.pop('--GCmodel', None)

        for option, argument in options.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                self.call.extend((option, argument))
        cmd = AtomicCmd(self.call, IN_BAM=options['BamPath'],
                        IN_BED=bed_path,
                        OUT_FILEPATH=destination+'.txt.gz')
        description = "<ANALYSIS:'%s', Bed:'%s', BAM: %s" % \
                      (anal_name, bed_name, base_bam)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class NucleoNode(CommandNode):
    def __init__(self, call, analysis_options, bedfile, dependencies=()):
        self.call = call
        bam_name, anal_name, outfmt, options = deepcopy(analysis_options)
        bed_name, bed_path = bedfile
        base_bam = os.path.split(bam_name)[-1]
        destination = \
            os.path.join(bam_name, outfmt(base_bam, bed_name, anal_name))
        if not options['Apply_GC_Correction']:
            options.pop('--GCmodel', None)

        for option, argument in options.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                self.call.extend((option, argument))
        cmd = AtomicCmd(self.call, IN_BAM=options['BamPath'],
                        IN_BED=bed_path,
                        OUT_FILEPATH=destination+'.txt.gz')
        description = "<ANALYSIS:'%s', Bed:'%s', BAM: %s" % \
                      (anal_name, bed_name, base_bam)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)


class GeneralExecuteNode(CommandNode):
    def __init__(self, call, analysis_options, bedfile, dependencies=()):
        self.call = call
        bam_name, anal_name, outfmt, options = deepcopy(analysis_options)
        bed_name, bed_path = bedfile
        base_bam = os.path.split(bam_name)[-1]
        destination = \
            os.path.join(bam_name, outfmt(base_bam, bed_name, anal_name))
        if not options.get('Apply_GC_Correction', True):
            options.pop('--GCmodel', None)


        for option, argument in options.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                self.call.extend((option, argument))
        cmd = AtomicCmd(self.call, IN_BAM=options['BamPath'],
                        IN_BED=bed_path,
                        OUT_FILEPATH=destination+'.txt.gz')
        description = "<ANALYSIS:'%s', Bed:'%s', BAM: %s" % \
                      (anal_name, bed_name, base_bam)
        CommandNode.__init__(self,
                             description=description,
                             command=cmd,
                             dependencies=dependencies)
