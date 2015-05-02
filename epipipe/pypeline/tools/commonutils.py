import pysam
import re
FMT = '{pref}{chr_no}'.format


def unpack(chrom, start, end, *rest):
    return chrom, start, end, rest


def corr_chrom(prefix, chrom):
    curr_chrom = ''.join(char for char in chrom if char.isdigit() or char in 'XYMT')
    return FMT(pref=prefix, chr_no=curr_chrom)


def read_bed(args):
    bampref = args.BamPrefix
    if args.bed:
        with open(args.bed, 'r') as bedfile:
            for line in bedfile:
                chrom, start, end, rest = unpack(*(re.split(r'\s+',
                                                 line.rstrip())))
                yield (corr_chrom(bampref, str(chrom)), int(start), int(end))


def read_mappa(args):
    bampref = args.BamPrefix
    with open(args.MappabilityPath, 'r') as mappafile:
        for line in mappafile:
            chrom, start, end, rest = unpack(*(re.split(r'\s+', line.rstrip())))
            yield (corr_chrom(bampref, str(chrom)),
                   int(start), int(end), float(rest[-1]))


class Cache(object):
    ''' class doc '''

    def __init__(self, filename, seq_len=1e6):
        self._fasta = pysam.Fastafile(filename)
        self._seq_len = int(seq_len)
        self._last_chrom = None
        self._fasta_str = None
        self._last_start = None
        self._actual_pos = None
        self._end = None

    def fetch_string(self, chrom, start, nbases):
        ''' docstring '''
        if self._last_chrom != chrom or (start-self._last_start) >= \
                self._seq_len or start >= self._end - nbases or \
                start < self._last_start:

            self._end = start + self._seq_len
            self._fasta_str = self._fasta.fetch(chrom,
                                                start=start, end=self._end)
            self._last_start = start
            self._last_chrom = chrom
        self._actual_pos = start-self._last_start
        return self._fasta_str[self._actual_pos:self._actual_pos+nbases]

    def closefile(self):
        ''' docstring '''
        return self._fasta.close()


class GC_correction(object):
    ''' class doc '''
    def _GCmodel_ini(self):
        if self.arg.GCmodel:
            self._fasta_dat = Cache(self.arg.FastaPath)
            with open(self.arg.GCmodel, 'r') as f:
                next(f)  # do not need the header
                self._model = [float(line.rstrip('\n').split('\t')[-1])
                               for line in f]
                self._GC_model_len = len(self._model)

    def _get_gc_corr_dep(self, pos):
        if self.arg.GCmodel:
            fasta_str = self._fasta_dat.fetch_string(self.chrom,
                                                     pos, self._GC_model_len-1)
            gc_idx = fasta_str.count('G')+fasta_str.count('C')
            return self._model[gc_idx]
        else:
            return 1
