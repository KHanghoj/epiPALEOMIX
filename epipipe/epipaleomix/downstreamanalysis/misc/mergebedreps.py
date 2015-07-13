from __future__ import print_function
import re, sys
FMT = '{c}\t{s}\t{e}\t{cov}\t{methyl}\n'.format

# takes only chromosomes without "chr"

def sum_dfs(df1, df2):
    methreads = (df1['cov']*float(df1['methyl']))+(df2['cov']*float(df2['methyl']))
    totreads = (df1['cov']+df2['cov'])
    df1['cov'] = totreads
    df1['methyl'] = methreads/totreads
    # df1['cov'] = (df1['cov']+df2['cov'])/2.0
    # df1['methyl'] = (df1['methyl']+df2['methyl'])/2.0
    return df1

def unpack(chrom, start, end, cov, methyl, *rest):
    return {'c':int(chrom), 's':int(start),'e':int(end),
            'cov':int(cov), 'methyl':float(methyl)}

def nextline(f):
    try:
        return unpack(*re.split(r'\s+', f.next().rstrip()))
    except StopIteration:
        return None

def finishedfile(fout, f):
    for line in f:
        fout.write(line)

outname, in1, in2 = sys.argv[1:]

print(outname, in1, in2,file=sys.stderr)
with open(outname, 'w') as fout:
    with open(in1, 'r') as fone:
        with open(in2, 'r') as ftwo:
            row_two = nextline(ftwo)
            row_one = nextline(fone)                        
            while True:
                if not row_one:                # finish row_two and quit
                    if row_two:
                        fout.write(FMT(**row_two))
                        finishedfile(fout, ftwo)
                    break

                if not row_two:                # finished row_one and quit
                    if row_one:
                        fout.write(FMT(**row_one))
                        finishedfile(fout, fone)
                    break

                if row_one['c'] == row_two['c']: ## chromosome
                    if row_one['s'] < row_two['s']:  # one smallest
                        fout.write(FMT(**row_one))
                        row_one = nextline(fone)                        
                    elif row_one['s'] > row_two['s']:  # two smallest
                        fout.write(FMT(**row_two))
                        row_two = nextline(ftwo)
                    else:
                        fout.write(FMT(**sum_dfs(row_one, row_two)))
                        row_two = nextline(ftwo)
                        row_one = nextline(fone)                        
                elif row_one['c'] > row_two['c']:
                    fout.write(FMT(**row_two))
                    row_two = nextline(ftwo)
                else:
                    fout.write(FMT(**row_one))
                    row_one = nextline(fone)  
