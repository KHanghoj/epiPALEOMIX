#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import time
import logging
import pypeline.yaml
import pypeline.logger
import optparse
from set_procname import set_procname
from nodes.gccorrect_Node import \
    GccorrectNode, CreateGCModelNode, \
    GccorrectNode_Mid  # , GccorrectNode_Final
from nodes.execute_Node import \
    GeneralExecuteNode, \
    General_Plot_Node
from nodes.cleanbedfiles_Node import CleanFilesNode
from pypeline.node import MetaNode
from epiomix_makefile import read_epiomix_makefile
from pypeline.pipeline import Pypeline
from pypeline.common.console import \
    print_err, \
    print_info

FINETUNERANGE = [-10, -5, 5, 10]


class MakefileError(RuntimeError):
    """Raised if a makefile is unreadable, or does not meet specifications."""


def parse_args(argv):
    ''' simpleparser '''
    usage_str = "%prog <command> [options] [makefiles]"
    # version_str = "%%prog %s" % (pypeline.__version__,)
    parser = optparse.OptionParser(usage=usage_str)
    # parser = optparse.OptionParser(usage=usage_str, version=version_str)
    parser.add_option('--temp-root', type=str, default='./temp')
    parser.add_option("--dry-run", action="store_true", default=False,
                      help="If passed, only a dry-run in performed"
                      ", the dependency ")
    parser.add_option('--destination', type=str, default='./OUTPUT')
    parser.add_option('--max-threads', type=int, default=4)
    pypeline.logger.add_optiongroup(parser)
    return parser.parse_args(argv)


def check_path(temp_dir):
    if not os.path.exists(temp_dir):
        try:
            os.makedirs(temp_dir)
        except OSError, error:
            print_err("ERROR: Could not create temp root:\n\t%s" % (error,))
            return 1


def filter_bedfiles(temp_root, d_make):
    bednodes = []
    if d_make.bedfiles.get('EnabledFilter', False):
        uniqueness = d_make.bedfiles.get('UniquenessFilter', 0)
        mappapath = d_make.prefix.get('--MappabilityPath', '')
        for bedname, in_bedpath in checkbedfiles_ext(d_make.bedfiles):
            out_bedpath = os.path.join(temp_root,
                                       os.path.basename(in_bedpath))
            d_make.bedfiles[bedname] = out_bedpath  # RENEW BEDFILE PATH
            bednodes.append(CleanFilesNode(mappapath, uniqueness,
                                           in_bedpath, out_bedpath))
    return bednodes


def split_bedfiles(temp_root, d_make, filternode=()):
    pass


def checkbedfiles_ext(bedfiles):
    for bedname, bedpath in bedfiles.iteritems():
        if isinstance(bedpath, str) and bedpath.endswith('.bed'):
            yield bedname, bedpath


def main_anal_to_run(opts):
    for analysis, options in opts.iteritems():
        if analysis in ANALYSES.keys() and options['Enabled']:
            yield analysis, ANALYSES[analysis]


def makegcnodes(d_bam, rang, subnodes=()):
    return [GccorrectNode(d_bam, rl, subnodes=subnodes) for rl in rang]


def calc_gcmodel(d_bam):
    if d_bam.opts['GCcorrect'].get('Enabled', False):
        nodes, lastnode = [], []
        scale = d_bam.opts['GCcorrect']['Resolution']
        rlmin, rlmax = \
            d_bam.opts['GCcorrect'].get('MapMinMaxReadLength', [56, 57])
        d_bam.opts['NucleoMap']['--MaxReadLen'] = rlmax
        nodes = makegcnodes(d_bam, xrange(rlmin, rlmax+1, scale))
        midnode = [GccorrectNode_Mid(d_bam, subnodes=nodes)]
        lastnode = makegcnodes(d_bam, FINETUNERANGE, subnodes=midnode)
        topnode = CreateGCModelNode(d_bam, subnodes=lastnode)
        gcoutfile = (out for out in topnode.output_files if
                     out.endswith('.txt')).next()
        d_bam.opts['BamInfo']['--GCmodel'] = gcoutfile
        return [topnode]
    return []


def run_analyses(anal, a_path, d_bam, d_make, bedinfo, m_node):
    aux_python, aux_R = a_path
    bed_name, bed_path = bedinfo
    node = GeneralExecuteNode(anal, aux_python, d_bam, bed_name, bed_path,
                              dependencies=m_node)
    # if bedplot is false -> no plot
    if aux_R == '_' or not d_make.bed_plot[bed_name]:
        return node
    infile = (out for out in node.output_files if out.endswith('txt.gz'))
    return General_Plot_Node(aux_R, infile.next(), anal, dependencies=node)


def make_metanode(depen_nodes, bamname):
    descrip_fmt = 'Metanode: {}: correction of GC and bedfile'.format
    return MetaNode(description=descrip_fmt(bamname),
                    dependencies=depen_nodes)


ANALYSES = {'Phasogram': ['./tools/phasogram.py', './tools/phaso.R'],
            'WriteDepth': ['./tools/pileupdepth.py', '_'],
            'NucleoMap': ['./tools/nucleomap.py', './tools/nucleo_merge.R'],
            'MethylMap': ['./tools/methylmap.py', './tools/methyl_merge.R']
            }


class makef_collect(object):
    def __init__(self, make):
        self.makefile = make.pop('Makefile', {})
        self.prefix = self.makefile.pop('Prefixes', {})
        self.beddata = self.makefile.pop('BedFiles', {})
        self.bedfiles, self.bed_plot = self._splitbedopts()

    def _splitbedopts(self):
        bedf, bedp = {}, {}
        for bedname, bedopts in self.beddata.iteritems():
            if isinstance(bedopts, dict):
                bedf[bedname] = bedopts["Path"]
                bedp[bedname] = bedopts["MakeMergePlot"]
            else:
                bedf[bedname] = bedopts
        return bedf, bedp


class bam_collect(object):
    def __init__(self, config, bam_name, opts, d_make):
        self.bam_name = bam_name
        self.i_path = os.path.join(config.temp_root, self.bam_name)
        self.o_path = os.path.join(config.destination, self.bam_name)
        self._createpaths()
        self.opts = opts
        self.baminfo = self.opts['BamInfo']
        self.fmt = '{}_{}_{}.txt.gz'
        self.prefix = d_make.prefix  # comes from makef_collect class
        self.generalopts = self._coerce_to_dic(self.baminfo, self.prefix)

    def _createpaths(self):
        check_path(self.i_path)
        check_path(self.o_path)

    def retrievedat(self, analtype):
        return self._coerce_to_dic(self.opts[analtype], self.generalopts)

    def _coerce_to_dic(self, *args):
        dic = {}
        for arg in args:
            if not isinstance(arg, dict):
                arg = dict([arg])
            dic.update(arg)
        return dic


def run(config, makefiles):
    check_path(config.temp_root)
    logfile_template = time.strftime("Epiomix_pipe.%Y%m%d_%H%M%S_%%02i.log")
    pypeline.logger.initialize(config, logfile_template)
    logger = logging.getLogger(__name__)
    pipeline = Pypeline(config)
    topnodes = []
    for make in read_epiomix_makefile(makefiles):
        d_make = makef_collect(make)
        filternode = filter_bedfiles(config, d_make)
        for bam_name, opts in d_make.makefile['BamInputs'].items():
            d_bam = bam_collect(config, bam_name, opts, d_make)
            gcnode = calc_gcmodel(d_bam)
            m_node = make_metanode(gcnode+filternode, d_bam.bam_name)
            for bedinfo in checkbedfiles_ext(d_make.bedfiles):
                for anal, a_path in main_anal_to_run(opts):
                    topnodes.append(run_analyses(anal, a_path, d_bam,
                                                 d_make, bedinfo, m_node))
            if not topnodes:
                topnodes.extend(gcnode+filternode)
    pipeline.add_nodes(topnodes)
    logger.info("Running BAM pipeline ...")
    pipeline.run(dry_run=config.dry_run, max_running=config.max_threads)


def _print_usage():
    basename = os.path.basename(sys.argv[0])
    if basename == "./epipaleomix.py":
        basename = "epipaleomix"

    print_info("Epipaleomix\n")
    print_info("Usage:")
    print_info("  -- %s help           -- Display this message" % basename)
    print_info("  -- %s dryrun [...]   -- Perform dry run of pipeline on provided makefiles." % basename)
    # print_info("     %s                   Equivalent to 'epipaleomix run --dry-run [...]'." % (" " * len(basename),))
    print_info("  -- %s run [...]      -- Run pipeline on provided makefiles." % basename)


def main(argv):
    set_procname("epipaleomix")
    try:
        config, args = parse_args(argv)
        if args and args[0].startswith("dry"):
            config.dry_run = True
    except RuntimeError, error:
        print_err(error)
        return 1

    commands = ("run", "dry_run", "dry-run", "dryrun")
    if (len(args) == 0) or (args[0] not in commands):
        _print_usage()
        return 1
    elif not args[1:]:
        _print_usage()
        print_err("\nPlease specify at least one makefile!")
        return 1
    # NOTE THAT WE SPLICE OUT ANOTHER INPUT
    # BEFORE PASSING ON TO RUN FUNC
    return run(config, args[1:])

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
