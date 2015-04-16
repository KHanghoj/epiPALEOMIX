from __future__ import print_function
import sys
import os
import time
import pypeline.yaml
import pypeline.logger
import optparse
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


def parse_args(argv):
    ''' simpleparser '''
    parser = optparse.OptionParser()
    parser.add_option('--temp-root', type=str, default='./temp')
    parser.add_option('--run', action='store_false', default=True)
    parser.add_option('--outputfolder', type=str, default='./OUTPUT')
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


def filter_bedfiles(bedfiles, destination_pref, mappapath):
    bednodes = []
    uniqueness = bedfiles.get('UniquenessFilter', 0)
    for bedname, in_bedpath in checkbedfiles_ext(bedfiles):
        out_bedpath = os.path.join(destination_pref,
                                   os.path.basename(in_bedpath))
        bedfiles[bedname] = out_bedpath  # renew the path to the bedfiles
        bednodes.append(CleanFilesNode(mappapath, uniqueness,
                                       in_bedpath, out_bedpath))
    return bednodes


def checkbedfiles_ext(bedfiles):
    for bedname, bedpath in bedfiles.iteritems():
        if isinstance(bedpath, str) and bedpath.endswith('.bed'):
            yield bedname, bedpath


def main_anal_to_run(opts):
    return {key: ANALYSES[key] for key, val in opts.iteritems() if
            key in ANALYSES.keys() and val['Enabled']}


def calc_gccorrectionmodel(opts, prefix_opt, io_paths):
    nodes = []
    GC_resolution = opts['GCcorrect']['Resolution']
    suffix = io_paths['bamname']+'_GCcorrect'
    gc_dic = coerce_to_dic(prefix_opt, opts['BamInfo'], opts['GCcorrect'])
    rlmin, rlmax = gc_dic.pop('MapMinMaxReadLength', [56, 57])

    for rl in xrange(rlmin, rlmax+1, GC_resolution):  # this is each readlength
        nodes.append(GccorrectNode(gc_dic, suffix, io_paths, rl))
    topnode = CreateGCModelNode(suffix, io_paths, subnodes=nodes)

    opts['NucleoMap']['--MaxReadLen'] = rlmax
    opts['BamInfo']['--GCmodel'] = list(topnode.output_files)[0]
    return [topnode]


def get_io_paths(config, bam_name):
    io_paths = {'o_out': os.path.join(config.outputfolder, bam_name),
                'temp': os.path.join(config.temp_root, bam_name),
                'bamname': bam_name}
    check_path(io_paths['o_out'])
    check_path(io_paths['temp'])
    return io_paths


def make_gcnode(opts, prefix_opt, io_paths):
    gc_bedfile_node = []
    if opts['GCcorrect'].get('Enabled', False):
        gc_bedfile_node.extend(calc_gccorrectionmodel(opts, prefix_opt,
                               io_paths))
    return gc_bedfile_node


def run_filterbed(bedfiles, config, prefix_opt):
    bedfilenode = []
    if bedfiles.get('EnabledFilter', False):
        bedfilenode.extend(filter_bedfiles(bedfiles,
                           config.temp_root,
                           prefix_opt["--MappabilityPath"]))
    return bedfilenode


def run_analyses(aux_paths, analysis_options, BedInfo, mNode):
    aux_python, aux_R = aux_paths
    node = GeneralExecuteNode(aux_python, analysis_options, BedInfo,
                              dependencies=mNode)
    curr_anal = analysis_options[1]
    if aux_R == '_':  # we have no plotting for writedepth
        return node
    return General_Plot_Node(aux_R, list(node.output_files)[0],
                             curr_anal, dependencies=node)


def makemetanode(gcnode, dependenNode, bamname):
    descrip_fmt = 'Metanode: {}: correction of GC and bedfile'.format
    mNode = MetaNode(description=descrip_fmt(bamname),
                     dependencies=dependenNode + gcnode)
    return mNode

ANALYSES = {'Phasogram': ['./tools/phasogram.py', './tools/phaso.R'],
            'WriteDepth': ['./tools/pileupdepth.py', '_'],
            'NucleoMap': ['./tools/nucleomap.py', './tools/nucleo_merge.R'],
            'MethylMap': ['./tools/methylmap.py', './tools/methyl_merge.R']
            }
GEN_OUTPUT_FMT = '{}_{}_{}'.format


def main(argv):
    config, makefiles = parse_args(argv)
    # this is not what we normally want:
    check_path(config.temp_root)
    logfile_template = time.strftime("Epiomix_pipe.%Y%m%d_%H%M%S_%%02i.log")

    pypeline.logger.initialize(config, logfile_template)
    pipeline = Pypeline(config)
    topnodes = []
    for make in read_epiomix_makefile(makefiles):
        dependenNode = []
        makefile = make.get('Makefile')
        prefix_opt = makefile.get('Prefixes')
        bedfiles = makefile.get('BedFiles')
        dependenNode.extend(run_filterbed(bedfiles, config, prefix_opt))
        for bam_name, opts in makefile['BamInputs'].items():
            io_paths = get_io_paths(config, bam_name)
            mNode = makemetanode(make_gcnode(opts, prefix_opt, io_paths),
                                 dependenNode, bam_name)
            bam_analyses = main_anal_to_run(opts)
            general_opts = coerce_to_dic(prefix_opt, opts["BamInfo"])
            for BedInfo in checkbedfiles_ext(bedfiles):
                for analysis, aux_paths in bam_analyses.iteritems():
                    anal_opt = [io_paths['o_out'], analysis, GEN_OUTPUT_FMT,
                                coerce_to_dic(general_opts, opts[analysis])]
                    topnodes.append(run_analyses(aux_paths, anal_opt,
                                                 BedInfo, mNode))
    pipeline.add_nodes(topnodes)
    pipeline.run(dry_run=config.run, max_running=config.max_threads)
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
