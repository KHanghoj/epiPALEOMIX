#!/usr/bin/env python
import os
import optparse
import pypeline.ui
from pypeline.config import PerHostConfig, PerHostValue, ConfigError

__version_info__ = (1, 0, 0)
__version__ = 'v%i.%i.%i' % __version_info__ 

__foreword__ = """
An automated pipeline for generating epigenomic 
profiles from high-throughput sequencing shotgun
 data of sequencing data underlying Ancient samples
===================================================
Kristian Hanghoej,(Jakob Skou etc.), Ludovic Orlando.
Centre for GeoGenetics
Natural History Museum of Denmark
University of Copenhagen.

Aarhus blabla

Version: %s
                              -----
""" %  __version__

__description__ = """%prog Automatic generation of phasogram,
nucleomsome and methylation maps using shotgun 
high through-put DNA sequencing underlying ancient samples provided in bam format.""" 

__epilog__ = """Please report bugs and suggestions for improvements to:
Kristian Hanghoej (kristianhanghoej@gmail.com)
Jakob and Tobias
Ludovic Orlando (lorlando@snm.ku.dk)
"""

class EpiPaleomixIndentHelper(optparse.IndentedHelpFormatter):
    def format_usage(self, usage):
        return usage + "\n"
    def format_epilog(self, epilog):
        return "\n" + epilog


def _run_config_parser(argv):
    # Helper class for optparse.OptionParser
    # allows reading from and writing to a config file
    per_host_cfg = PerHostConfig("epiPALEOMIX")
    usage_str = "%prog <command> [options] [makefiles]"
    version_str  = "%%prog %s" % (__version__,)
    parser = optparse.OptionParser(usage = __foreword__ + "\n" + usage_str,
                                   version = version_str,
                                   description = __description__,
                                   epilog = __epilog__,
                                   formatter=EpiPaleomixIndentHelper(),
                                   prog="epiPALEOMIX")

    pypeline.ui.add_optiongroup(parser,
                                ui_default=PerHostValue("quiet"),
                                color_default=PerHostValue("on"))
    pypeline.logger.add_optiongroup(parser, default = PerHostValue("warning"))

    parser.add_option("--max-threads", type = int, default = per_host_cfg.max_threads,
                     help = "Maximum number of threads to use in total [%default]")
    parser.add_option("--dry-run", action = "store_true", default = False,
                     help = "If passed, only a dry-run in performed, the dependency "
                     "tree is printed, and no tasks are executed.")
    parser.add_option("--temp-root", default = per_host_cfg.temp_root,
                     help = "Location for temporary files and folders [%default/]")
    parser.add_option("--destination", default = './',
                     help = "The destination folder for result files. By default, "
                     "files will be placed in ./OUT_{makefile name}/")
    return per_host_cfg.parse_args(parser, argv)


def parse_config(argv):
    config, args = _run_config_parser(argv)
    pypeline.ui.set_ui_colors(config.ui_colors)
    return config, args
