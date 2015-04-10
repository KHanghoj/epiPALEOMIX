import sys
import os
import time
import pypeline.yaml
import pypeline.logger
import optparse
from gccorrect_Node import \
    GccorrectNode, \
    CreateGCModelNode
from pypeline.common.console import print_err
from epiomix_makefile import read_epiomix_makefile
from pypeline.pipeline import Pypeline


def parse_args(argv):
    parser = optparse.OptionParser()
    parser.add_option('--temp-root', default=os.path.join(os.getcwd(), 'temp'))
    pypeline.logger.add_optiongroup(parser)

    return parser.parse_args(argv)


class analyses_dic_arguments(object):
    """docstring for analyses_dic_arguments"""
    def __init__(self, arg):
        super(analyses_dic_arguments, self).__init__()
        self.arg = arg


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


def main_good(argv):
    # lst = []
    config, makefiles = parse_args(argv)
    # this is not what we normally want:
    check_path(config.temp_root)
    nodes = []
    logfile_template = time.strftime("Epiomix_pipeline."
                                     "%Y%m%d_%H%M%S_%%02i.log")

    pypeline.logger.initialize(config, logfile_template)
    pipeline = Pypeline(config)
    for make in read_epiomix_makefile(makefiles):
        # lst.append(copy.deepcopy(make))
        makefile = make.pop('Makefile', {})
        # stats = make.pop('Statistics', {})

        fasta_opt = makefile.pop('Prefixes', {})
        for each_bam, opts in makefile['BamInputs'].iteritems():
            gc_dic = {}
            if opts['GCcorrect'].pop('Enabled', False):
                updatedic(gc_dic, fasta_opt,
                          opts['BamInfo'],
                          opts['GCcorrect'])
                nodes.append(GccorrectNode(gc_dic, config.temp_root))
            print(nodes)
    pipeline.add_nodes(nodes)
    pipeline.run(dry_run=False, max_running=2)


def dest_prefix(*args):
    return(os.path.join(*args))


def main(argv):
    # lst = []
    config, makefiles = parse_args(argv)
    # this is not what we normally want:
    check_path(config.temp_root)

    logfile_template = time.strftime("Epiomix_pipeline."
                                     "%Y%m%d_%H%M%S_%%02i.log")

    pypeline.logger.initialize(config, logfile_template)
    pipeline = Pypeline(config)
    for make in read_epiomix_makefile(makefiles):
        # lst.append(copy.deepcopy(make))
        makefile = make.pop('Makefile', {})
        # stats = make.pop('Statistics', {})
        topnodes = []
        fasta_opt = makefile.pop('Prefixes', {})
        bednodes = func()
        for each_bam, opts in makefile['BamInputs'].iteritems():
            nodes = []
            bamfile = opts['BamInfo'].pop("BamPath")
            b_name = os.path.basename(bamfile)
            gc_dic = {}
            if opts['GCcorrect'].pop('Enabled', False):
                upd_dic(gc_dic, fasta_opt, opts['BamInfo'], opts['GCcorrect'])

                destination_pref = dest_prefix(config.temp_root, b_name+'_')

                rlmin, rlmax = gc_dic.pop('MapMinMaxReadLength', [56, 57])
                for rl in xrange(rlmin, rlmax+1, 1):  # this is each readlength
                    nodes.append(GccorrectNode(gc_dic, bamfile,
                                               destination_pref, rl))
                topnodes.append(CreateGCModelNode(config.temp_root, b_name,
                                os.getcwd(), dependencies=nodes))
            # filternodes = list()
            # for bedfile, path in beds.iteritems():
            #     newpath = os.path.basename(path)
            #     shutil.copy(path, os.path.join(config.temp_root,newpath)
            #     update makefile[bedfile] = newpath
            #     filternodes.append(FILTERNODE(Rscript)))
            # loop all bedfiles. copy to new temp folder. replace
            # path ned to new one and filter with mappability script
            # then make a meta node of the two top nodes.
            # simple
            # dernaest skal vi starte de 4/5 individuelle analyser
            # dog skal vi igennem bedfile to gange. dog ikke 
            # saa vigtigt tror jeg.

    pipeline.add_nodes(topnodes)
    pipeline.run(dry_run=False, max_running=2)
if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
