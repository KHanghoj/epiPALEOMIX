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
    parser.add_option('--dry-run', action='store_true', default=False)
    parser.add_option('--outputfolder', type=str, default='./OUTPUT')
    parser.add_option('--max-threads', type=int, default=2)
    pypeline.logger.add_optiongroup(parser)
    return parser.parse_args(argv)


def check_path(temp_dir):
    if not os.path.exists(temp_dir):
        try:
            os.makedirs(temp_dir)
        except OSError, error:
            print_err("ERROR: Could not create temp root:\n\t%s" % (error,))
            return 1


def upd_dic(dic, *args):
    for arg in args:
        dic.update(arg)


def coerce_to_dic(*args):
    dic = {}
    for arg in args:
        if not isinstance(arg, dict):
            arg = dict([arg])
        dic.update(arg)
    return dic


# def dest_prefix(*args):
#     return(os.path.join(*args))


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


def analyses_to_run(opts):
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


def run_gc(opts, prefix_opt, io_paths):
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

# def nucleofunc(NucleoNode):
    #
    # return RUNRSCRIPTNODE(dependencies=NucleoNode)
    # append this function to topnodes and done deal
    # within the loop.
NucleoFunc, MethylFunc, PhasoFunc, WriteFunc =\
    object(), object(), object(), object()

# ANALYSES = {'NucleoMap': NucleoFunc,
#             'MethylMap': MethylFunc,
#             'Phasogram': PhasoFunc,
#             'WriteDepth': WriteFunc}


def main_anal(aux_paths, analysis_options, BedInfo, mNode):
    aux_python, aux_R = aux_paths
    node = GeneralExecuteNode(aux_python, analysis_options, BedInfo,
                              dependencies=mNode)
    infile = list(node.output_files)[0]
    if aux_R == '_':
        return node
    return General_Plot_Node(aux_R, infile, analysis_options[1],
                             dependencies=node)


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
            descrip_fmt = 'Metanode: {}: correction of GC and bedfile'.format
            mNode = MetaNode(description=descrip_fmt(io_paths['bamname']),
                             dependencies=dependenNode +
                             run_gc(opts, prefix_opt, io_paths))
            bam_analyses = analyses_to_run(opts)
            general_opts = coerce_to_dic(prefix_opt, opts["BamInfo"])
            # for BedInfo in checkbedfiles_ext(bedfiles):
            #     for analysis, aux_path in bam_analyses.iteritems():
            #         analysis_options = \
            #             [io_paths['o_out'], analysis, GEN_OUTPUT_FMT,
            #              coerce_to_dic(general_opts, opts[analysis])]
            #         topnodes.append(GeneralExecuteNode(aux_path,
            #                         analysis_options, BedInfo,
            #                         dependencies=mNode))
            for BedInfo in checkbedfiles_ext(bedfiles):
                for analysis, aux_paths in bam_analyses.iteritems():
                    analysis_options = \
                        [io_paths['o_out'], analysis, GEN_OUTPUT_FMT,
                         coerce_to_dic(general_opts, opts[analysis])]
                    topnodes.append(main_anal(aux_paths, analysis_options,
                                              BedInfo, mNode))
    pipeline.add_nodes(topnodes)
    pipeline.run(dry_run=config.dry_run, max_running=config.max_threads)
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
