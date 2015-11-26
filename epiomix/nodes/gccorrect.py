#!/usr/bin/env python
from __future__ import print_function
from pypeline.node import CommandNode, Node
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.common.fileutils import move_file, reroot_path
import pypeline.common.versions as versions

from epiomix.tools import gccorrect_mid, gccorrect

import os
Rscript_VERSION = versions.Requirement(call   = ("Rscript", "--version"),
                                       search = r"version (\d+)\.(\d+)\.(\d+)",
                                       checks = versions.GE(2, 15, 3))

GC_NAME = '_GCcorrect'
### this is for individual readlength gccorrection
class GccorrectNode(Node):
    def __init__(self, d_bam, rl, dependencies=()):
        self.dest = os.path.join(d_bam.bam_temp_local,
                                 d_bam.bam_name+GC_NAME+'_'+str(rl))
        self.rl, self.d_bam  = rl, d_bam
        description = ("<Gccorrect: '%s' window length: '%s' based on chromosome '%s' >" %
                       (self.dest, rl, self.d_bam.opts['GCcorrect']['--ChromUsed']))

        Node.__init__(self,
                      description=description,
                      input_files=self.d_bam.baminfo["BamPath"],
                      output_files=self.dest,
                      dependencies=dependencies)
        assert len(self.output_files) == 1, self.output_files

    def _run(self, _config, temp):
        dest = reroot_path(temp, self.dest)
        self.inputs = [self.d_bam.baminfo["BamPath"], dest]
        self._add_options('GCcorrect')
        self.inputs.extend(["--ReadLength", str(self.rl)])
        gccorrect.main(self.inputs)

    def _teardown(self, _config, temp):
        move_file(reroot_path(temp, self.dest), self.dest)
        Node._teardown(self, _config, temp)

    def _add_options(self, name):
        optargs = self.d_bam.retrievedat(name)
        for option, argument in optargs.iteritems():
            if isinstance(option, str) and option.startswith('-'):
                if not isinstance(argument, str):
                    argument = str(argument)
                self.inputs.extend((option, argument))


class CreateGCModelNode(CommandNode):
    def __init__(self, d_bam, dependencies=()):
        aux_r = os.path.join(os.path.dirname(gccorrect.__file__),
                             'model_gc.R')
        
        call = ['Rscript', aux_r, '%(IN_SOURCE)s', str(d_bam.bam_name+GC_NAME),
                str(len(dependencies)), '%(OUT_FILEPATH)s', '%(OUT_PLOT)s']
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
                             dependencies=dependencies)
