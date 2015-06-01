#!/usr/bin/env python
from __future__ import print_function
import sys
import os
import time
import logging
import pypeline.yaml
import pypeline.logger
import optparse
from epipaleomix.epi_mkfile import epi_mkfile
from epipaleomix.tools import checkchromprefix
from epipaleomix.set_procname import set_procname
from epipaleomix.epi_mkfile.epi_makefile import read_epiomix_makefile
from epipaleomix.nodes.gccorrect_Node import \
    GccorrectNode, \
    CreateGCModelNode, \
    GccorrectNode_Mid
from epipaleomix.nodes.execute_Node import \
    GeneralExecuteNode, \
    General_Plot_Node
from epipaleomix.nodes.cleanbedfiles_Node import \
    CleanFilesNode, \
    SplitBedFile, \
    MergeDataFiles
from pypeline.node import MetaNode
from pypeline.pipeline import Pypeline
from pypeline.common.console import \
    print_err, \
    print_info


FINETUNERANGE = [-10, -5, 5, 10]
ANALYSES = ['Phasogram', 'WriteDepth', 'NucleoMap', 'MethylMap']


class MakefileError(RuntimeError):
    """Raised if a makefile is unreadable, or does not meet specifications."""


class make_collect(object):
    def __init__(self, make):
        self.makefile = make.pop('Makefile', {})
        self.prefix = self.makefile.pop('Prefixes', {})
        self.beddata = self.makefile.pop('BedFiles', {})
        self.bedfiles, self.bed_plot = self._splitbedopts()

    def _splitbedopts(self):
        ''' Split bedfiles options between bedname path and plot boolean '''
        bedf, bedp = {}, {}
        for bedn, bedopts in self.beddata.iteritems():
            if isinstance(bedopts, dict):
                bedf[bedn] = bedopts["Path"]
                bedp[bedn] = bedopts["MakeMergePlot"]
            else:  # just options to the specific bedfile
                bedf[bedn] = bedopts
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
        self.prefix = d_make.prefix  # comes from make_collect class

    def _createpaths(self):
        check_path(self.i_path)
        check_path(self.o_path)

    def retrievedat(self, anal):
        return self._coerce_to_dic(self.opts[anal], self.baminfo,
                                   self.prefix)

    def _coerce_to_dic(self, *args):
        dic = {}
        for arg in args:
            if not isinstance(arg, dict):
                arg = dict([arg])
            dic.update(arg)
        return dic


def parse_args(argv):
    ''' simpleparser '''
    usage_str = "%prog <command> [options] [makefiles]"
    version_str = "%%prog %s" % (pypeline.__version__,)
    parser = optparse.OptionParser(usage=usage_str, version=version_str)
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


def split_bedfiles(config, d_make):
    uniqueness = d_make.bedfiles.get('UniquenessFilter', 0)
    mappapath = d_make.prefix.get('--MappabilityPath', False)
    enabl_filter = d_make.bedfiles.get('EnabledFilter', False)
    filtnode, nodes = [], []
    for bedn, in_bedp in checkbedfiles_ext(d_make.bedfiles):
        if enabl_filter and mappapath:
            filtnode = [CleanFilesNode(config, in_bedp, mappapath, uniqueness)]
            d_make.bedfiles[bedn] = ''.join(filtnode[0].output_files)
        splnode = SplitBedFile(config, d_make.bedfiles, bedn, subnodes=filtnode)
        nodes.append(splnode)
    return nodes


def checkbedfiles_ext(bedfiles):
    for bedname, bedpath in bedfiles.items():
        if isinstance(bedpath, str) and bedpath.endswith('.bed'):
            yield bedname, bedpath


def checkbed_list(bedfiles):
    for bedname, bedpaths in bedfiles.items():
        if (isinstance(bedpaths, list) and
            all([bedp.endswith('.bed') for bedp in bedpaths])):
            yield bedname, bedpaths


def main_anal_to_run(opts):
    for analysis, options in opts.iteritems():
        if analysis in ANALYSES and options['Enabled']:
            yield analysis


def makegcnodes(d_bam, rang, subnodes=()):
    return [GccorrectNode(d_bam, rl, subnodes=subnodes) for rl in rang]


def concat_gcsubnodes(nodecls, bam, ran, subn=()):
    return [nodecls(bam, subnodes=makegcnodes(bam, ran, subn))]


def calc_gcmodel(d_bam):
    if d_bam.opts['GCcorrect'].get('Enabled', False):
        rlmin, rlmax = \
            d_bam.opts['GCcorrect'].get('MapMinMaxReadLength', MakefileError)
        return concat_gcsubnodes(CreateGCModelNode, d_bam, FINETUNERANGE,
                                   subnodes=concat_gcsubnodes(GccorrectNode_Mid,
                                                              d_bam, xrange(rlmin, rlmax+1, 15)))
    return []


def run_analyses(anal, d_bam, d_make, bedinfo, m_node):
    bedn, bed_paths = bedinfo
    nodes = []
    for idx, bed_p in enumerate(bed_paths):
        nodes.append(GeneralExecuteNode(anal, d_bam, bedn+str(idx), bed_p, dependencies=m_node))
    mergenode = MergeDataFiles(d_bam, anal, bedn, subnodes=nodes)
    if not d_make.bed_plot[bedn] or anal == 'WriteDepth': # if bedplot is false -> no plot
        return mergenode
    infile = (out for out in mergenode.output_files).next()
    return General_Plot_Node(infile, anal, dependencies=[mergenode])


def make_metanode(depen_nodes, bamname):
    descrip_fmt = "Metanode: '{}': GC correction and bedfile split"
    return MetaNode(description=descrip_fmt.format(bamname),
                    dependencies=depen_nodes)


def check_chrom_prefix(bedfiles, d_make):
    for bam_name, opts in d_make.makefile['BamInputs'].items():
        baminfo = opts['BamInfo']
        for bedfile, bedpath in checkbedfiles_ext(bedfiles):
            checkchromprefix.main([baminfo['BamPath'],
                                   d_make.prefix.get('--FastaPath'),
                                   bedpath])


def run(config, makefiles):
    check_path(config.temp_root)
    logfile_template = time.strftime("Epipaleomix_pipe.%Y%m%d_%H%M%S_%%02i.log")
    pypeline.logger.initialize(config, logfile_template)
    logger = logging.getLogger(__name__)
    pipeline = Pypeline(config)
    topnodes = []
    for make in read_epiomix_makefile(makefiles):
        d_make = make_collect(make)
        check_chrom_prefix(d_make.bedfiles, d_make)
        splitbednode = split_bedfiles(config, d_make)
        for bam_name, opts in d_make.makefile['BamInputs'].items():
            d_bam = bam_collect(config, bam_name, opts, d_make)
            gcnode = calc_gcmodel(d_bam)
            m_node = make_metanode(gcnode+splitbednode, d_bam.bam_name)
            for bedinfo in checkbed_list(d_make.bedfiles):
                for anal in main_anal_to_run(opts):
                    topnodes.append(run_analyses(anal, d_bam,
                                                 d_make, bedinfo, m_node))
            if not topnodes:
                topnodes.extend(gcnode+splitbednode)
    pipeline.add_nodes(topnodes)
    logger.info("Running Epipaleomix pipeline ...")
    pipeline.run(dry_run=config.dry_run, max_running=config.max_threads)


def _print_usage():
    basename = os.path.basename(sys.argv[0])
    if basename == "./epipaleomix.py":
        basename = "epipaleomix"
    print_info("Epipaleomix\n")
    print_info("Usage:")
    print_info("  -- %s help           -- Display this message" % basename)
    print_info("  -- %s dryrun [...]   -- Perform dry run of pipeline on provided makefiles." % basename)
    print_info("  -- %s makefile [...] -- Generate makefile template to std.out." % basename)
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
    commands = ("makefile", "mkfile", "run", "dry_run", "dry-run", "dryrun")
    if (len(args) == 0) or (args[0] not in commands):
        _print_usage()
        return 1
    elif args[0] in ("mkfile", "makefile"):
        return epi_mkfile.main(args[1:])
    elif not args[1:]:
        _print_usage()
        print_err("\nPlease specify at least one makefile!")
        return 1
    # NOTE THAT WE SPLICE OUT ANOTHER INPUT BEFORE PASSING ON TO RUN FUNC
    return run(config, args[1:])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
