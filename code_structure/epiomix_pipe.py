from __future__ import print_function
# NucleoNode, MethylNode, PhasoNode, WriteNode = object(), object(), object(), object()
import sys
import os
import time
import pypeline.yaml
import pypeline.logger
import optparse
from nodes.gccorrect_Node import \
    GccorrectNode, CreateGCModelNode
from nodes.execute_Node import \
    PhasoNode, \
    WriteNode, \
    NucleoNode
from nodes.cleanbedfiles_Node import CleanFilesNode
from pypeline.node import MetaNode
from pypeline.common.console import print_err
from epiomix_makefile import read_epiomix_makefile
from pypeline.pipeline import Pypeline


def parse_args(argv):
    parser = optparse.OptionParser()
    parser.add_option('--temp-root', type=str, default='./temp')
    parser.add_option('--outputfolder', type=str, default='./output')
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
    uniqueness = bedfiles.pop('UniquenessFilter')
    # uniqueness = bedfiles.get('UniquenessFilter', None)
    for bedname, in_bedpath in checkbedfiles_ext(bedfiles):
        out_bedpath = os.path.join(destination_pref,
                                   os.path.basename(in_bedpath))
        bedfiles[bedname] = out_bedpath  # renew the path to the bedfiles
        bednodes.append(CleanFilesNode(mappapath, uniqueness,
                                       in_bedpath, out_bedpath))
    return bednodes


def checkbedfiles_ext(bedfiles):
    for bedname, bedpath in bedfiles.items():
        if isinstance(bedpath, str) and bedpath.endswith('.bed'):
            yield bedname, bedpath


def analyses_to_run(opts):
    return {key: ANALYSES[key] for key, val in opts.iteritems() if
            key in ANALYSES.keys() and val['Enabled']}


def calc_gccorrectionmodel(opts, fasta_opt, io_paths):
    nodes, gc_dic = [], {}
    suffix = io_paths['bamname']+'_GCcorrect'
    temp_dest_pref = os.path.join(io_paths['temp'], suffix)
    gc_dic = coerce_to_dic(fasta_opt, opts['BamInfo'], opts['GCcorrect'])
    rlmin, rlmax = gc_dic.pop('MapMinMaxReadLength', [56, 57])
    for rl in xrange(rlmin, rlmax+1, 5):  # this is each readlength
        nodes.append(GccorrectNode(gc_dic, temp_dest_pref, rl))
    opts['BamInfo']['--GCmodel'] = \
        os.path.join(io_paths['o_out'], 'GC_Model_%s.txt' % (str(suffix)))
    return [CreateGCModelNode(suffix, io_paths,
            dependencies=nodes)]


# def nucleofunc(NucleoNode):
    #
    # return RUNRSCRIPTNODE(dependencies=NucleoNode)
    # append this function to topnodes and done deal
    # within the loop.



# ANALYSES = {'NucleoMap': NucleoFunc,
#             'MethylMap': MethylFunc,
#             'Phasogram': PhasoFunc,
#             'WriteDepth': WriteFunc}

# ANALYSES = {'NucleoMap': NucleoNode,
#             'MethylMap': MethylNode,
#             'Phasogram': PhasoNode,
#             'WriteDepth': WriteNode}

ANALYSES = {'Phasogram': PhasoNode,
            'WriteDepth': WriteNode,
            'NucleoMap': NucleoNode}

EXECUTABLES = {
    'NucleoMap': ['python', './tools/nucleomap.py', '%(IN_BAM)s',
                  '%(IN_BED)s', '%(OUT_FILEPATH)s'],
    'MethylMap': ['python', './tools/methylmap.py', '%(IN_BAM)s',
                  '%(IN_BED)s', '%(OUT_FILEPATH)s'],
    'Phasogram': ['python', './tools/phasogram.py', '%(IN_BAM)s',
                  '%(IN_BED)s', '%(OUT_FILEPATH)s'],
    'WriteDepth': ['python', './tools/pileupdepth.py', '%(IN_BAM)s',
                   '%(IN_BED)s', '%(OUT_FILEPATH)s']
}


def get_io_paths(config, bam_name):
    io_paths = {'o_out': os.path.join(config.outputfolder, bam_name),
                'temp': os.path.join(config.temp_root, bam_name),
                'bamname': bam_name}
    check_path(io_paths['o_out'])
    check_path(io_paths['temp'])
    return io_paths


def main(argv):
    config, makefiles = parse_args(argv)
    # this is not what we normally want:
    check_path(config.temp_root)
    logfile_template = time.strftime("Epiomix_pipe.%Y%m%d_%H%M%S_%%02i.log")

    pypeline.logger.initialize(config, logfile_template)
    pipeline = Pypeline(config)
    topnodes = []
    for make in read_epiomix_makefile(makefiles):
        makefile = make.pop('Makefile')
        fasta_opt = makefile.pop('Prefixes')
        bedfilenode = []
        if makefile['BedFiles'].get('EnabledFilter', False):
            bedfilenode.extend(filter_bedfiles(makefile['BedFiles'],
                               config.temp_root,
                               fasta_opt["--MappabilityPath"]))

        for bam_name, opts in makefile['BamInputs'].items():
            gc_bedfile_node = []
            io_paths = get_io_paths(config, bam_name)
            if opts['GCcorrect'].get('Enabled', False):
                gc_bedfile_node.extend(calc_gccorrectionmodel(opts, fasta_opt,
                                       io_paths))
            mNode = MetaNode(description="Metanode: '{}': for correction of "
                                         "GC and bedfiles".format
                                         (io_paths['bamname']),
                             dependencies=gc_bedfile_node+bedfilenode)
            bam_analyses = analyses_to_run(opts)
            outputfmt = '{}_{}_{}'.format
            general_opts = coerce_to_dic(fasta_opt, opts["BamInfo"])
            for BedInfo in checkbedfiles_ext(makefile['BedFiles']):
                for analysis, Nodeclass in bam_analyses.iteritems():
                    call = EXECUTABLES[analysis]
                    analysis_options = [io_paths['o_out'], analysis, outputfmt,
                                        coerce_to_dic(general_opts,
                                                      opts[analysis])]
                    topnodes.append(Nodeclass(call, analysis_options,
                                    BedInfo, dependencies=mNode))
    pipeline.add_nodes(topnodes)
    pipeline.run(dry_run=True, max_running=config.max_threads)
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
