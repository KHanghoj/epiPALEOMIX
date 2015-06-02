from __future__ import print_function
import gzip, re, sys, argparse
phaserange = 250
FMT='{}\t{}\t{}\t{}\t{}\t{}\n'.format
def unpackepipal(chrom, start, end, depth, score, bedc):
    return chrom, int(start), int(end), float(depth), float(score), str(bedc)


def p_a(args):
    parser = argparse.ArgumentParser(prog='chunk nucleosome data')
    parser.add_argument('nuclpath', type=str)
    parser.add_argument('--maxdyaddist', default=250, type=int)
    parser.add_argument('--mindyaddist', default=140, type=int)
    parser.add_argument('--mincoverage', default=1, type=int)
    parser.add_argument('--scorecutoff', default=1, type=int)
    return parser.parse_args(args)


def writetofile(Dat):
    chrom, last_start, runningsum = Dat.out()
    sys.stdout.write(FMT(chrom, last_start[0], last_start[-1],
                         len(last_start), runningsum,
                         runningsum/len(last_start)))


class datacontain(object):
    def __init__(self):
        self.last_score = self.runscore = 0
        self.last_chrom = ''
        self.last_start = []

    def reset(self, chrom, start, score):
        self.last_score = self.runscore = score
        self.last_chrom = chrom
        self.last_start = [start]

    def scoreupdate(self, start, score, former=None):
        if former:
            self.runscore += score-self.last_score
            self.last_start[-1] = start
        else:
            self.runscore += score
            self.last_start.append(start)
        self.last_score = score

    def out(self):
        return self.last_chrom, self.last_start, self.runscore

        
def nuclchunks(args):
    MAXDYADDIST = args.maxdyaddist
    MINDYADDIST = args.mindyaddist
    COVERAGE = args.mincoverage
    SCORE = args.scorecutoff
    with gzip.open(args.nuclpath) as f_in:
        f_in.next()  # removes header
        Dat = datacontain()
        for line in f_in:
            chrom, start, end, depth, score, bedc = unpackepipal(*re.split(r'\s+', line.rstrip()))
            if score>=SCORE:                
                if chrom != Dat.last_chrom:
                    if len(Dat.last_start) >= COVERAGE:
                        writetofile(Dat)
                    Dat.reset(chrom, start, score)

                if start-Dat.last_start[-1] <= MAXDYADDIST:
                    if start-Dat.last_start[-1] >= MINDYADDIST:
                        Dat.scoreupdate(start, score)
                    elif Dat.last_score < score:
                        Dat.scoreupdate(start, score, True)
                else:
                    if len(Dat.last_start) >= COVERAGE:
                        writetofile(Dat)
                    Dat.reset(chrom, start, score)
        if len(Dat.last_start) >= COVERAGE:
            writetofile(Dat)


def main(argv):
    args = p_a(argv)
    nuclchunks(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


