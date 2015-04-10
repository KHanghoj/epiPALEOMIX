from pypeline.node import CommandNode
from pypeline.atomiccmd.command import AtomicCmd
from pypeline.atomiccmd.sets import ParallelCmds
import tempfile
import os
import shutil


class CleanFilesNode(CommandNode):
    def __init__(self, makefile, destfile, dependencies=()):
        
        pass

        pass
        # CommandNode.__init__(self,
        #                      description=description,
        #                      command=cmd,
        #                      dependencies=dependencies)
