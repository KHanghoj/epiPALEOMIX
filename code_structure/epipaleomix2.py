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
    GccorrectNode, CreateGCModelNode
from nodes.execute_Node import \
    GeneralExecuteNode, \
    General_Plot_Node
from nodes.cleanbedfiles_Node import CleanFilesNode
from pypeline.node import MetaNode
from pypeline.common.console import print_err
from epiomix_makefile import read_epiomix_makefile
from pypeline.pipeline import Pypeline
from pypeline.common.console import \
    print_err, \
    print_info


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


def coerce_to_dic(*args):
    dic = {}
    for arg in args:
        if not isinstance(arg, dict):
            arg = dict([arg])
        dic.update(arg)
    return dic


def filter_bedfiles(temp_root, dat_make):
    bednodes = []
    uniqueness = dat_make.bedfiles.get('UniquenessFilter', 0)
    mappapath = dat_make.prefix.get('--MappabilityPath', '')
    for bedname, in_bedpath in checkbedfiles_ext(dat_make.bedfiles):
        out_bedpath = os.path.join(temp_root,
                                   os.path.basename(in_bedpath))
        dat_make.bedfiles[bedname] = out_bedpath  # RENEW BEDFILE PATH
        bednodes.append(CleanFilesNode(mappapath, uniqueness,
                                       in_bedpath, out_bedpath))
    return bednodes


def checkbedfiles_ext(bedfiles):
    for bedname, bedpath in bedfiles.iteritems():
        if isinstance(bedpath, str) and bedpath.endswith('.bed'):
            yield bedname, bedpath


# def main_anal_to_run(opts):
#     return {key: ANALYSES[key] for key, val in opts.iteritems() if
#             key in ANALYSES.keys() and val['Enabled']}


def main_anal_to_run(opts):
    for analysis, options in opts.iteritems():
        if analysis in ANALYSES.keys() and options['Enabled']:
            yield analysis, ANALYSES[analysis]


def calc_gcmodel(dat_bam, dat_make):
    nodes = []
    scale = dat_bam.opts['GCcorrect']['Resolution']
    rlmin, rlmax = \
        dat_bam.opts['GCcorrect'].get('MapMinMaxReadLength', [56, 57])

    for rl in xrange(rlmin, rlmax+1, scale):  # this is each readlength
        nodes.append(GccorrectNode(dat_bam, dat_make, rl))
    topnode = CreateGCModelNode(dat_bam, subnodes=nodes)

    dat_bam.opts['NucleoMap']['--MaxReadLen'] = rlmax
    gcoutfile = (out for out in topnode.output_files if out.endswith('.txt'))
    dat_bam.opts['BamInfo']['--GCmodel'] = gcoutfile.next()
    return [topnode]


def get_io_paths(config, bam_name):
    io_paths = {'o_out': os.path.join(config.destination, bam_name),
                'temp': os.path.join(config.temp_root, bam_name),
                'bamname': bam_name}
    check_path(io_paths['o_out'])
    check_path(io_paths['temp'])
    return io_paths


def make_gcnode(dat_bam, dat_make):
    if dat_bam.opts['GCcorrect'].get('Enabled', False):
        return calc_gcmodel(dat_bam, dat_make)
    return []


def run_filterbed(config, dat_make):
    if dat_make.bedfiles.get('EnabledFilter', False):
        # return filter_bedfiles(bedfiles, config.temp_root,
        #                        prefix["--MappabilityPath"])
        return filter_bedfiles(config.temp_root, dat_make)
    return []


def run_analyses(anal, a_path, dat_bam, bedinfo, m_node):
    aux_python, aux_R = a_path
    node = GeneralExecuteNode(anal, aux_python, dat_bam, bedinfo,
                              dependencies=m_node)
    if aux_R == '_':  # we have no plotting for writedepth
        return node
    infile = (out for out in node.output_files if
              out.endswith('txt.gz')).next()
    return General_Plot_Node(aux_R, infile, anal, dependencies=node)


def make_metanode(depen_nodes, bamname):
    descrip_fmt = 'Metanode: {}: correction of GC and bedfile'.format
    return MetaNode(description=descrip_fmt(bamname),
                    dependencies=depen_nodes)


ANALYSES = {'Phasogram': ['./tools/phasogram.py', './tools/phaso.R'],
            'WriteDepth': ['./tools/pileupdepth.py', '_'],
            'NucleoMap': ['./tools/nucleomap.py', './tools/nucleo_merge.R'],
            'MethylMap': ['./tools/methylmap.py', './tools/methyl_merge.R']
            }
GEN_OUTPUT_FMT = '{}_{}_{}'.format


class makef_collect(object):
    def __init__(self, make):
        self.makefile = make.pop('Makefile', {})
        self.prefix = self.makefile.pop('Prefixes', {})
        self.bedfiles = self.makefile.pop('BedFiles', {})


class bam_collect(object):
    def __init__(self, config, bam_name, opts, dat_make):
        self.bam_name = bam_name
        self.i_path = os.path.join(config.temp_root, self.bam_name)
        self.o_path = os.path.join(config.destination, self.bam_name)
        self.baminfo = opts['BamInfo']
        self.opts = opts
        self.dat_make = dat_make  # comes from makef_collect class
        self.fmt = GEN_OUTPUT_FMT

    def retrievedat(self, analtype):
        # potentially retrieve data from makef_collect prefix.
        return self.opts[analtype]+self.baminfo+self.dat_make.prefix


def run(config, makefiles):
    # config, makefiles = parse_args(argv)
    # this is not what we normally want:
    check_path(config.temp_root)
    logfile_template = time.strftime("Epiomix_pipe.%Y%m%d_%H%M%S_%%02i.log")
    pypeline.logger.initialize(config, logfile_template)
    logger = logging.getLogger(__name__)
    pipeline = Pypeline(config)
    topnodes = []
    for make in read_epiomix_makefile(makefiles):
        dat_make = makef_collect(make)
        # makefile = make.get('Makefile')
        # prefix = makefile.get('Prefixes')
        # bedfiles = makefile.get('BedFiles')
        # filternode.extend(run_filterbed(bedfiles, config, prefix))
        filternode = run_filterbed(config, dat_make)
        for bam_name, opts in dat_make.makefile['BamInputs'].items():
            dat_bam = bam_collect(config, bam_name, opts, dat_make)
            # io_paths = get_io_paths(config, bam_name)
            m_node = make_metanode(make_gcnode(dat_bam, dat_make) +
                                   filternode)
            # general_opts = coerce_to_dic(prefix, opts["BamInfo"])
            for bedinfo in checkbedfiles_ext(dat_make.bedfiles):
                for anal, a_path in main_anal_to_run(opts):
                    # anal_opt = [io_paths['o_out'], anal, GEN_OUTPUT_FMT,
                    #             coerce_to_dic(general_opts, opts[anal])]
                    t_node = run_analyses(anal, a_path, dat_bam, bedinfo, m_node)
                    topnodes.append(run_analyses(a_path, anal_opt,
                                                 bedinfo, m_node))
    pipeline.add_nodes(topnodes)
    logger.info("Running BAM pipeline ...")
    pipeline.run(dry_run=config.dry_run, max_running=config.max_threads)


def _print_usage():
    basename = os.path.basename(sys.argv[0])
    if basename == "./epipaleomix.py":
        basename = "epipaleomix"

    # print_info("Epipaleomix %s\n" % (pypeline.__version__,))
    print_info("Epipaleomix\n")
    print_info("Usage:")
    print_info("  -- %s help           -- Display this message" % basename)
    # print_info("  -- %s makefile [...] -- Generate makefile from 'SampleSheet.csv' files." % basename)
    print_info("  -- %s dryrun [...]   -- Perform dry run of pipeline on provided makefiles." % basename)
    print_info("     %s                   Equivalent to 'epipaleomix run --dry-run [...]'." % (" " * len(basename),))
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
    # Note that we remove antother input before passing on to run func
    return run(config, args[1:])

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
