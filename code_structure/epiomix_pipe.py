from __future__ import print_function
import sys
import os
import time
import pypeline.yaml
import pypeline.logger
import optparse
from nodes.gccorrect_Node import \
    GccorrectNode, CreateGCModelNode
from nodes.CleanBedfilesNode import CleanFilesNode
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


def dest_prefix(*args):
    return(os.path.join(*args))


def filter_bedfiles(bedfiles, destination_pref, mappapath):
    bednodes = []
    uniqueness = bedfiles.pop('UniquenessFilter')
    for bedname, in_bedpath in bedfiles.iteritems():
        if in_bedpath.endswith('.bed'):
            out_bedpath = os.path.join(destination_pref,
                                       os.path.basename(in_bedpath))
            bedfiles[bedname] = out_bedpath  # renew the path to the bedfiles
            bednodes.append(CleanFilesNode(mappapath, uniqueness,
                                           in_bedpath, out_bedpath))
    return bednodes


def calc_gccorrectionmodel(bamfile, opts, fasta_opt, temp_root, bamname):
    nodes, gc_dic = [], {}
    temp_gcpath = dest_prefix(temp_root, bamname+'_GCcorrect')
    gc_dic = coerce_to_dic(fasta_opt, opts['BamInfo'], opts['GCcorrect'])
    rlmin, rlmax = gc_dic.pop('MapMinMaxReadLength', [56, 57])
    for rl in xrange(rlmin, rlmax+1, 5):  # this is each readlength

        nodes.append(GccorrectNode(gc_dic, bamfile,
                                   temp_gcpath, rl))
    opts['BamInfo']['--GCmodel'] = os.path.join(bamname,
                                                'GC_Model_%s.txt'
                                                % (str(bamname)))
    return [CreateGCModelNode(temp_gcpath, bamname,
            dependencies=nodes)]

ANALYSES = ('NucleoMap', 'MethylMap', 'Phasogram', 'WriteDepth')
EXECUABLE_DIC = {
    'NucleoMap': ['python', '.tools/nucleomap.py', '%(IN_BAM)s',
                  '%(OUT_MERGE)s', '%(OUT_RAW)s'],
    'MethylMap': ['python', '.tools/methylmap.py', '%(IN_BAM)s',
                  '%(OUT_MERGE)s', '%(OUT_RAW)s'],
    'Phasogram': ['python', '.tools/phasogram.py', '%(IN_BAM)s',
                  '%(OUT_STD)s'],
    'WriteDepth': ['python', '.tools/pileupdepth.py', '%(IN_BAM)s',
                   '%(OUT_STD)s']
}


def main(argv):
    config, makefiles = parse_args(argv)
    # this is not what we normally want:
    check_path(config.temp_root)
    # check_path(config.outputfolder)
    logfile_template = time.strftime("Epiomix_pipe.%Y%m%d_%H%M%S_%%02i.log")

    pypeline.logger.initialize(config, logfile_template)
    pipeline = Pypeline(config)
    topnodes = []
    for make in read_epiomix_makefile(makefiles):
        makefile = make.pop('Makefile')
        fasta_opt = makefile.pop('Prefixes')
        bednodes = []
        if makefile['BedFiles'].pop('EnabledFilter', False):
            bednodes = filter_bedfiles(makefile['BedFiles'],
                                       config.temp_root,
                                       fasta_opt["--MappabilityPath"])
        else:  # REMOVE UNIQUENESS FILTER TO AVOID CONFLICT IF NOT USED
            makefile['BedFiles'].pop('UniquenessFilter')

        for bam_name, opts in makefile['BamInputs'].items():
            check_path(bam_name)

            bamfile = opts['BamInfo'].pop('BamPath')
            if opts['GCcorrect'].pop('Enabled', False):
                gcnodes = calc_gccorrectionmodel(bamfile, opts, fasta_opt,
                                                 config.temp_root, bam_name)
                topnodes = gcnodes+bednodes if bednodes else gcnodes
                # GCmodel is in the makefile now
            mNode = MetaNode(description="'%s': metanode for correction of "
                                         "GC and bedfiles" % (bam_name, ),
                             dependencies=topnodes)
            print(mNode)

            # # put this mnode into every script from now on.
            # # make a meta node here of all the files above and start
            # # the actual analysis.
            # # dependenciNode = MetaNode()
            bam_analyses = [key for key, val in opts.iteritems() if
                            key in ANALYSES and
                            val['Enabled']]
            for bedfile in makefile['BedFiles'].items():
                dic = coerce_to_dic(fasta_opt, bedfile, opts["BamInfo"])
                for analysis in bam_analyses:
                    # input:
                    dic = EXECUABLE_DIC[analysis]
                    print(dic)
                    # bamfile, bedfile,outputpaths
                    # mNode to each analysis

    pipeline.add_nodes(topnodes)
    pipeline.run(dry_run=False, max_running=config.max_threads)
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
