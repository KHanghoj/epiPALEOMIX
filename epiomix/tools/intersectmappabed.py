from __future__ import print_function
import argparse
import re
import sys
import os

BEDFMT = '{}_{}_{}_{}'.format
FMT = '{}\t{}\t{}\t{}\n'.format


def parse_args(argv):
    ''' docstring '''
    parser = argparse.ArgumentParser(prog='GCcorrection')
    parser.add_argument('bed', type=str)
    parser.add_argument('MappabilityPath', type=str)
    parser.add_argument('MappaUniqueness', help="default is 0.9", type=float)
    # parser.add_argument('out', nargs='?', default=sys.stdout)
    parser.add_argument('out')
    return parser.parse_known_args(argv)


def getstrand(rest):
    curr_strand = \
        re.search(r"(\s+|^)(?P<strandtype>[+-])(\s+|$)", " ".join(rest))
    if curr_strand:
        return curr_strand.group("strandtype")
    return '+'


def getbed(args):
    checkfile(args.bed)
    with open(args.bed, 'r') as data:
        d = {}
        for line in data:
            c, s, e, rest = unpack(*(re.split(r'\s+', line.rstrip())))
            bedcoord = (BEDFMT(c, s, e, getstrand(rest)), )
            if c in d.keys():
                d[c].append((int(s), int(e), bedcoord))
            else:
                d[c] = [(int(s), int(e), bedcoord)]

    for chrom in d.keys():
        d[chrom].sort()
        # tuplesorted = sorted(d[chrom])
        # d[chrom] = tuplesorted

    return d


def unpack(chrom, start, end, *rest):
    return chrom, int(start), int(end), rest


def _read_mappa(args):
    checkfile(args.MappabilityPath)
    with open(args.MappabilityPath, 'r') as mappafile:
        for line in mappafile:
            chrom, start, end, rest = unpack(*(re.split(r'\s+',
                                                        line.rstrip())))
            if float(rest[-1]) >= args.MappaUniqueness:
                yield str(chrom), int(start), int(end)


def getmappa(args):
    d = {}
    l = []
    cprev = ''
    endprev = -1
    for chrom, start, end in _read_mappa(args):

        if chrom != cprev and l:
            endprev = end
            d[cprev] = sorted(l)
            l = []

        if start <= endprev and l:
            l[-1][-1] = end  # update previous end
            cprev = chrom
            endprev = end
            continue

        l.append([start, end])
        cprev = chrom
        endprev = end

    if l:
        d[cprev] = sorted(l)
    return d


def setupcallregion(chrom, handler):
    def output(regbegin, regend, rest):
        handler.write(FMT(chrom, regbegin, regend, "\t".join(map(str, rest))))
    return output


def getmappacoord(mappaiter):
    try:
        return mappaiter.next()
    except StopIteration:
        # if outside of mappabilty list, high values of
        # start and end are assigned
        # to run out of bed file. Might be error prone
        return int(1e10), int(1e11)


def getintersect(fout, bed, mappa, chrom):
    callregion = setupcallregion(chrom, fout)
    bedtarget = bed[chrom]

    mappaiter = iter(mappa[chrom])
    mappa_start, mappa_end = getmappacoord(mappaiter)
    lastend = -1
    for start, end, rest in bedtarget:
        # restart mappa if new bed region overlaps previous
        if start <= lastend:
            mappaiter = iter(mappa[chrom])
            mappa_start, mappa_end = getmappacoord(mappaiter)

        # bed end is prior to mappability start, new bedfiles
        if end < mappa_start:
            continue
        # bed start adhead of mappability end, new mapparegion
        while start > mappa_end:
            mappa_start, mappa_end = getmappacoord(mappaiter)

        while end > mappa_end:  # print out regions
            callregion(max(start, mappa_start), min(end, mappa_end), rest)
            mappa_start, mappa_end = getmappacoord(mappaiter)
        if end > mappa_start:
            callregion(max(start, mappa_start), min(end, mappa_end), rest)
        lastend = end


def checkfile(x):
    assert os.path.exists(x), "file does not exists: '%s'" % (x,)


def getoverlaps(args):
    bed_regions = getbed(args)
    mappability_regions = getmappa(args)
    assert all((True if k in mappability_regions else False
                for k in bed_regions)), ("Not all chromosomes in the bedfile "
                                         "are not present in Mappability file")
    with open(args.out, 'w') as fout:
        for chrom in sorted(bed_regions):
            getintersect(fout, bed_regions, mappability_regions, chrom)


def main(argv):
    args, unknown = parse_args(argv)
    getoverlaps(args)
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
