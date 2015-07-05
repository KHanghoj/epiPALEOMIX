from __future__ import print_function
import gzip, re, sys
FMT = '{c}\t{s}\t{d}\t{t}\t{bedc}\n'.format
### THIS merging file takes a out and two methylation files as input. It is expected that the files are numerically sorted in following order: bedc 1_10000_11000, chrom, start
### bedc must be splitted first to 1 10000 11000 then sort

def sum_dfs(df1, df2):
    df1['d'] = df1['d']+df2['d']
    df1['t'] = df1['t']+df2['t']
    return df1

def unpack(chrom, start, dea, tot, bedc):
    return {'c':str(chrom), 's':int(start),
            'd':int(dea), 't':int(tot), 'bedc':str(bedc)}

def nextline(f):
    try:
        return unpack(*re.split(r'\s+', f.next().rstrip()))
    except StopIteration:
        return None

def bedcdict(chrom, start, end):
    return {'c':int(chrom.replace('chr','')), 's':int(start), 'e':int(end)}
    
    
def checkbedc(fout, row_one, row_two):
    flag=0
    one = bedcdict(*row_one['bedc'].split('_'))
    two = bedcdict(*row_two['bedc'].split('_'))
    if one['c'] < two['c']:  # chrom one smallest
        fout.write(FMT(**row_one))
        flag=1
    elif one['c'] > two['c']:  # two smallest
        fout.write(FMT(**row_two))
        flag=2
    else:
        if one['s'] < two['s']:
            fout.write(FMT(**row_one))
            flag=1
        else:
            fout.write(FMT(**row_two))
            flag=2
    return flag

def finishedfile(fout, f):
    for line in f:
        fout.write(line)

outname, in1, in2 = sys.argv[1:]

print(outname, in1, in2,file=sys.stderr)
with gzip.open(outname, 'wb') as fout:
    with gzip.open(in1, 'rb') as fone:
        fout.write(fone.next()) # adding header to new merge fil
        ##fone.next() # remove header
        with gzip.open(in2, 'rb') as ftwo:
            ftwo.next() # remove header
            row_two = nextline(ftwo)
            row_one = nextline(fone)                        

            while True:

                if not row_one:                # finished row_two and quit
                    fout.write(FMT(**row_two))
                    finishedfile(fout, ftwo)
                    break

                if not row_two:                # finished row_one and quit
                    fout.write(FMT(**row_one))
                    finishedfile(fout, fone)
                    break

                if row_one['bedc'] == row_two['bedc']:
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
                else:
                    flag = checkbedc(fout, row_one,row_two)   # write the row with 'smallest region'
                    if flag==1:
                        row_one = nextline(fone)
                    elif flag==2:
                        row_two = nextline(ftwo)
                    else:
                        sys.exit('major error')
