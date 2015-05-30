from __future__ import print_function
import gzip, re, sys, argparse
phaserange = 250
FMT='{}\t{}\t{}\t{}\t{}\t{}\n'.format
def unpackepipal(chrom, start, end, depth, score, bedc):
    return chrom, int(start), int(end), float(depth), float(score), str(bedc)


def p_a(args):
    parser = argparse.ArgumentParser(prog='chunk nucleosome data')
    parser.add_argument('nuclpath', type=str)
    parser.add_argument('--phaserange', default=250, type=int)
    parser.add_argument('--mindyaddist', default=140, type=int)
    parser.add_argument('--mincoverage', default=1, type=int)
    parser.add_argument('--scorecutoff', default=1, type=int)
    return parser.parse_args(args)


def writetofile(chrom, last_start, runningsum):
    sys.stdout.write(FMT(chrom, last_start[0], last_start[-1], len(last_start), runningsum,  runningsum/len(last_start)))


class datacontain(object):
    def __init__(self, chrom, start, score):
        self.last_score = self.runscore = 0
        self.last_chrom = ''
        self.last_start = []

    def reset(self, chrom, start, score):
        self.last_score = self.runscore = score
        self.last_chrom = chrom
        self.last_start = [start]

    def scoreupdate(self, start, score, former=None):

        if former:
            self.runscore += score-former
            self.last_start[-1] = start
        else:
            self.runscore += score
            self.last_start.append(start)
        self.last_score = score

        
def nuclchunks(args):
    PHASERANGE = args.phaserange
    MINDYADDIST = args.mindyaddist
    COVERAGE = args.mincoverage
    SCORE = args.scorecutoff
    with gzip.open(args.nuclpath) as f_in:
        f_in.next()
        last_score, runningscore = 0, 0
        chrom_last, last_start = '', []
        Dat = datacontain()
        for line in f_in:
            chrom, start, end, depth, score, bedc = unpackepipal(*re.split(r'\s+', line.rstrip()))
            if score>=SCORE:                
                if chrom != Dat.last_chrom:
                    if len(last_start) >= COVERAGE:
                        writetofile(chrom_last, last_start, runningscore)
                    Dat.reset(chrom, start, score)

                if start-Dat.last_start[-1] <= PHASERANGE:
                    if start-Dat.last_start[-1] >= MINDYADDIST:
                        Dat.scoreupdate(start, score)
                        # last_start.append(start)
                        # Dat.runscore += score
                        # Dat.last_score = score
                    elif last_score < score:
                        Dat.scoreupdate(start, score, last_score)
                        # last_start[-1] = start
                        # runningscore += (score-last_score)
                        # last_score = score
                        # got this far
                else:
                    if len(last_start) >= COVERAGE:
                        writetofile(chrom, last_start, runningscore)
                    last_start = [start]
                    last_score = runningscore = score
        if len(last_start) >= COVERAGE:
            writetofile(chrom, last_start, runningscore)


def main(argv):
    args = p_a(argv)
    nuclchunks(args)
#    run = methylepipal if args.analysis == 'epipaleomix' else RRBS_anal
#    run(args)
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


