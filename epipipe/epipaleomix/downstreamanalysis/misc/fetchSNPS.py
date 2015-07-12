from __future__ import print_function
def getchrom(f):
    '''
    test = getchrom('snppos.txt')
    test.send(None)
    chrom = test.send(2)  # now chrom is all snps on chromosome 2
    ## never call same chromosome twice. 
    '''
    d=set()
    chrom = 0 
    with open(f, 'r') as fsnp:
        # some initialization
        while True:
            try:
                
                currchrom, coord = (int(n) for n in next(fsnp).split('_'))
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
