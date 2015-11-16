#!/usr/bin/env python
from __future__ import print_function
import operator
import re
import sys
import os
import copy
import time
import logging
import pypeline.yaml
import pypeline.logger
import optparse
from pypeline.node import MetaNode
from pypeline.pipeline import Pypeline
from pypeline.common.console import \
    print_err, \
    print_info
from epipaleomix.tools.commonutils import check_path
from epipaleomix.config import parse_config, __version__
from epipaleomix.epimakefile import epicreatemkfile
from epipaleomix.epimakefile.epivalidmkfile import read_epiomix_makefile
from epipaleomix.tools import checkchromprefix, \
    checkmappabilitychrom, \
    getminmax

from epipaleomix.tools.bamdatastructure import BamCollect, \
    MakeCollect, \
    MakefileError
from epipaleomix.nodes.gccorrect import \
    GccorrectNode, \
    CreateGCModelNode, \
    GccorrectMidNode
from epipaleomix.nodes.execute import \
    GeneralExecuteNode, \
    GeneralPlotNode
from epipaleomix.nodes.cleanbedfiles import \
    CleanFilesNode, \
    SplitBedFileNode, \
    MergeDataFilesNode

FINETUNERANGE = [-10, -5, 5, 10]
ANALYSES = ['Phasogram', 'WriteDepth', 'NucleoMap', 'MethylMap']


def check_bed_exist(config, infile):
    ''' Check if bedfiles have been split already, then use them
        As a changes in no. of threads from run to run is error prone
        in this setup'''
    filena, fileext = os.path.splitext(os.path.basename(infile))
    reg = r'{}_([0-9]+).bed'.format(filena)
    bedfiles=((f, int(re.search(reg,f).groups()[0])) for f in os.listdir(config.temp_local) if re.search(reg,f))
    ##  this is just cool dense code::
    ##  map(lambda x:x[0],sorted(a,key=itemgetter(1)))
    ##  return map(operator.itemgetter(0),sorted(bedfiles,key=operator.itemgetter(1)))
    return [os.path.join(config.temp_local, path) for path, val in sorted(bedfiles,key=operator.itemgetter(1))]


def split_bedfiles(config, d_make):
    uniqueness = d_make.bedfiles.get('UniquenessFilter', 0)
    mappapath = d_make.prefix.get('--MappabilityPath', False)
    enabl_filter = d_make.bedfiles.get('EnabledFilter', False)
    filtnode, nodes = [], []
    for bedn, in_bedp in checkbedfiles_ext(d_make.bedfiles):
        if enabl_filter and mappapath:
            filtnode = [CleanFilesNode(config, d_make, bedn, mappapath, uniqueness)]
        bedexists = check_bed_exist(config, d_make.bedfiles[bedn])
        if bedexists:
            d_make.bedfiles[bedn] = bedexists
        else:
            splnode = SplitBedFileNode(config, d_make, bedn, subnodes=filtnode)
            nodes.append(splnode)
    return nodes


def chromused_coerce_to_string(bam):
    chrused = bam.opts['GCcorrect'].get('--ChromUsed', MakefileError)
    bam.opts['GCcorrect']['--ChromUsed'] =  str(chrused)
    noregions = bam.opts['GCcorrect'].get('--NoRegions', MakefileError)
    if isinstance(noregions, str) and noregions.lower() == 'all':
        noregions = int(1e7)
    elif isinstance(noregions, int):
        noregions = int(noregions)
    else:
        raise MakefileError('--NoRegions: "%s" is incorrect. Must be of a string of {All, all, ALL} or an positive integer' % (noregions,))
    bam.opts['GCcorrect']['--NoRegions'] = noregions
    

def checkbedfiles_ext(bedfiles):
    for bedname, bedpath in bedfiles.items():
        if isinstance(bedpath, str) and bedpath.endswith('.bed'):
            yield bedname, bedpath


def checkbed_list(bedfiles):
    for bedname, bedpaths in bedfiles.items():
        if (isinstance(bedpaths, list) and
            all([bedp.endswith('.bed') for bedp in bedpaths])):
            yield bedname, bedpaths


def main_anal_to_run(bedinfo , opts):
    bedn, bed_paths = bedinfo    
    for analysis, options in opts.iteritems():
        if analysis in ANALYSES and options['Enabled'] and (bedn not in options['ExcludeBed']):
            yield analysis

def update_excludebed(d_make, d_bam):
    fmt = "{}".format
    if d_make.bedfiles.get('EnabledFilter', False):
        fmt = "{}MappaOnly".format

    for anal, opts in d_bam.opts.iteritems():
        if anal in ANALYSES:
            excl_bed = opts.get('ExcludeBed')
            if isinstance(excl_bed, str):
                opts['ExcludeBed'] = [fmt(excl_bed)]
            elif isinstance(excl_bed, list):
                opts['ExcludeBed'] = [fmt(bed) for bed in excl_bed]
            elif excl_bed is None:
                opts['ExcludeBed'] = []
            else:
                raise MakefileError('Exclude bed in %s is incorrect. Must be a str, list, or None' % (anal,))


def getdequelen(d_bam):
    rlmin, rlmax, gcmax = getminmax.main(d_bam.baminfo)
    d_bam.opts['WriteDepth']['--DequeLength'] = rlmax
    d_bam.opts['NucleoMap']['--DequeLength'] = rlmax
    return rlmin, gcmax

def concat_gcsubnodes(nodecls, d_bam, gcwindows, subn=()):
    return [nodecls(d_bam, subnodes=[GccorrectNode(d_bam, rl, subnodes=subn) for rl in gcwindows])]

### if individual readlength to be used
### changes nodes/gccorrect to used the nodes in the very bottom.
### used model_gc_individualreadlength.R instead of model_gc.R in createmodelNode 
### inherit indinode instead of nornormal in tools.commonutils

def calc_gcmodel(d_bam):
    rlmin, rlmax = getdequelen(d_bam)
    if d_bam.opts['GCcorrect'].get('Enabled', False):
        chromused_coerce_to_string(d_bam)
        checkmappabilitychrom.main([d_bam.prefix.get('--MappabilityPath', MakefileError),
                                    d_bam.opts['GCcorrect'].get('--ChromUsed', MakefileError)])
        resolution = 5
        return concat_gcsubnodes(CreateGCModelNode, d_bam,
                                 xrange(rlmin, rlmax+resolution, resolution))
    return []
    #     return concat_gcsubnodes(CreateGCModelNode,
    #                              d_bam,
    #                              FINETUNERANGE,
    #                              subn=concat_gcsubnodes(GccorrectMidNode, d_bam,
    #                                                     xrange(rlmin, rlmax+1, 15)))
    # return []

def check_chrom_prefix(d_make):
    for bam_name, opts in d_make.makefile['BamInputs'].items():
        baminfo = opts['BamInfo']
        for bedfile, bedpath in checkbedfiles_ext(d_make.bedfiles):
            checkchromprefix.main([baminfo['BamPath'],
                                   d_make.prefix.get('--FastaPath'),
                                   bedpath])


def run_analyses(anal, d_bam, d_make, bedinfo, splitbednode, gcnode):
    bedn, bed_paths = bedinfo
    nodes = []
    for idx, bed_p in enumerate(bed_paths):
        bedn_temp = '{}_0{}'.format(bedn, idx)
        nodes.append(GeneralExecuteNode(anal, d_bam, bedn_temp, bed_p,
                                        subnodes=gcnode,
                                        dependencies=splitbednode))
    mergenode = MergeDataFilesNode(d_bam, anal, bedn, subnodes=nodes)
    ## TODO:::: if writedepth and bedoptionmerge == TRUE: 
    ## TODO::::     apply the merger script that i just already on the mergeNode
    ## TODO:::: if methylation cleaning True in bed file. remove all SNPs based on databased already created
    ## TODO:::: if nucleomap and Advancednucleomap is TRUE. apply it here. 
    if not d_make.bed_plot[bedn] or anal == 'WriteDepth': # if bedplot is false -> no plot
        return mergenode
    infile = (out for out in mergenode.output_files).next()
    return GeneralPlot(infile, anal, dependencies=[mergenode])

def make_outputnames(config, make):
    filename = make["Statistics"]["Filename"]
    outputname = os.path.splitext(os.path.basename(filename))[0]
    config.makefiledest = os.path.join(config.destination, 'OUT_' + outputname)
    # config.temp_local = os.path.join(config.makefiledest, 'TEMPORARYFILES_' + outputname)
    config.temp_local = os.path.join(config.destination, 'TEMPORARYFILES_' + outputname)    
    check_path(config.makefiledest)
    check_path(config.temp_local)

    
def create_nodes(config, makefiles):
    topnodes = []
    for makefile in read_epiomix_makefile(makefiles):
        make_outputnames(config, makefile)
        d_make = MakeCollect(makefile)
        check_chrom_prefix(d_make)
        splitbednode = split_bedfiles(config, d_make)
        for bam_name, opts in d_make.makefile['BamInputs'].items():
            d_bam = BamCollect(config, bam_name, opts, d_make)
            gcnode = calc_gcmodel(d_bam)
            update_excludebed(d_make, d_bam)
            for bedinfo in checkbed_list(d_make.bedfiles):
                for anal in main_anal_to_run(bedinfo, opts):
                    topnodes.append(run_analyses(anal, d_bam,
                                                 d_make, bedinfo,
                                                 splitbednode, gcnode))
            if not topnodes:
                topnodes.extend(splitbednode+gcnode)
    return topnodes
    
def run(config, makefiles):
    check_path(config.temp_root)
    logfile_template = time.strftime("Epipaleomix_pipe.%Y%m%d_%H%M%S_%%02i.log")
    pypeline.logger.initialize(config, logfile_template)
    logger = logging.getLogger(__name__)
    pipeline = Pypeline(config=config)
    topnodes = create_nodes(config, makefiles)

    assert topnodes, "Pipeline found no analyses to run. Check %s" % (makefiles)
    pipeline.add_nodes(topnodes)
    if config.list_output_files:
        logger.info("Printing output files ...")
        pipeline.print_output_files()
        return 0
    elif config.list_executables:
        logger.info("Printing required executables ...")
        pipeline.print_required_executables()
        return 0

    logger.info("Running Epipaleomix pipeline ...")
    if not pipeline.run(dry_run=config.dry_run,
                        max_running=config.max_threads,
                        progress_ui=config.progress_ui):
        return 1
    return 0


def _print_usage():
    basename = os.path.basename(sys.argv[0])
    if basename == "./epipaleomix.py":
        basename = "epipaleomix"
    print_info("epiPALEOMIX Pipeline %s\n" %(__version__,))
    print_info("epiPALEOMIX\n")
    print_info("Usage:")
    print_info("  -- %s help           -- Display this message" % basename)
    print_info("  -- %s dryrun [...]   -- Perform dry run of pipeline on provided makefiles." % basename)
    print_info("  -- %s makefile [...] -- Generate makefile template to std.out." % basename)
    print_info("  -- %s run [...]      -- Run pipeline on provided makefiles." % basename)


def main(argv):
    try:
        config, args = parse_config(argv)
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
        return epicreatemkfile.main(args[1:])
    elif not args[1:]:
        _print_usage()
        print_err("\nPlease specify at least one makefile!")
        return 1
    return run(config, args[1:])


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
