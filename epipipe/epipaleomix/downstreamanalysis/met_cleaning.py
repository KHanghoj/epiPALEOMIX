from __future__ import print_function
import gzip, re, sys, argparse
FMT='{}\t{}\t{}\t{}\t{}\n'.format

def getSNPchrom(f):
    '''
    test = getchrom('snppos.txt')
    test.send(None)
    chrom = test.send(2)  # now chrom is all snps on chromosome 2
    ## never call same chromosome twice. 
    FORMAT OF SNP DATA Frame: 1_100200
    '''
    d=set()
    chrom = 0 
    with gzip.open(f, 'rb') as fsnp:
        # some initialization
        while True:
            try:
                currchrom, coord = map(int, fsnp.next().split('\t',2)[:2])
                ## currchrom, coord = (int(n) for n in fsnp.next().split('\t',2)[:2])
                ## currchrom, coord = (int(n) for n in next(fsnp).split('_'))
            except StopIteration:
                if d:
                    yield d
                    d.clear()
                else:
                    yield None
                break
            
            if currchrom == chrom:
                #d.add((chrom, coord)) # this for testing
                d.add(coord)
            elif currchrom > chrom:
                ## print('before: %s' % chrom)
                chrom = yield d
                # after a 'send' command it runs through the function until just
                # after "d" (so it return d) on this line and stops until next call that
                # reassigns chrom and proceses until it comes back to this line where it
                # return d and waits for the next "send"
                ## print('after: %s' % chrom)
                d.clear()
                lastchrom, lastcoord = currchrom, coord
                if chrom == lastchrom:
                    #d.add((lastchrom, lastcoord))
                    d.add(lastcoord)


def unpackepipal(chrom, start, dea, tot, bedc):
    return chrom, int(start), int(dea), int(tot), str(bedc)

def p_a(args):
    parser = argparse.ArgumentParser(prog='clean methylation data')
    parser.add_argument('metpath')
    parser.add_argument('out')
    return parser.parse_args(args)

def methylepipal(args):
    with gzip.open(args.metpath, 'rb') as f_in:
        f_in.next() # removes the header
        ds, cs, tot = [], [], []
        checked = 0
        ## fetchSNP = getSNPchrom('/home/krishang/data/SNP142/snppos.txt')
        fetchSNP = getSNPchrom('/home/krishang/data/SNP142/cpgrelatedSNP_allmut.txt.gz')

        fetchSNP.send(None)  # initialize
        h = ('#chrom\tpos\tdeaminatedsites\tcoverage\tbedcoord')

        with gzip.open(args.out, 'wb') as f_out:
            f_out.write(FMT(*h.split('\t')))
            lastchrom = ''
            for line in f_in:
                if checked % 500000 == 0:
                    print('checked {} methylated sites'.format(checked), file=sys.stderr)
                chrom, start, dea, TplusC, bedc = unpackepipal(*re.split(r'\s+', line.rstrip()))
                if lastchrom != chrom:
                    if chrom.replace('chr','').upper() == 'MT':
                        SNPset=set()
                    else:
                        SNPset = fetchSNP.send(int(chrom.replace('chr','')))   # works only on autosome chromo
                    lastchrom = chrom
                if start in SNPset or (dea > 3 and (float(dea)/TplusC) >= 0.33): # remove potential SNP's C>T. Does not represent methylatio
                    continue
                f_out.write(FMT(chrom,start,dea, TplusC,bedc))
                checked += 1
                lastchrom = chrom
                

def main(argv):
    args = p_a(argv)
    methylepipal(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

