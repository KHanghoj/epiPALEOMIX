import pysam


def read_bed_W(args):
    if args.bed:
        with open(args.bed, 'r') as bedfile:
            for line in bedfile:
                input_line = (line.rstrip('\n')).split('\t')[:3]
                chrom, start, end = input_line
                yield (str(chrom), int(start), int(end))


def read_bed_WO(args):
    if args.bed:
        with open(args.bed, 'r') as bedfile:
            for line in bedfile:
                input_line = (line.rstrip('\n')).split('\t')[:3]
                chrom, start, end = input_line
                yield (str(chrom.replace('chr', '')), int(start), int(end))


def strtobool(val):
    if isinstance(val, bool) or isinstance(val, int):
        return(val)
    val = val.lower()
    if val in ('y', 'yes', 't', 'true', 'on', '1'):
        return 1
    elif val in ('n', 'no', 'f', 'false', 'off', '0'):
        return 0
    else:
        raise ValueError("invalid truth value %r" % (val,))


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
    """docstring for GC_correction"""
    def _GCmodel_ini(self):
        if self.arg.GCmodel:
            with open(self.arg.GCmodel, 'r') as f:
                self._model = [float(line.rstrip('\n').split('\t')[-1])
                               for line in f]
                self._GC_model_len = len(self._model)

    def _get_gc_corr_dep(self, pos):
        if self.arg.GCmodel:
            fasta_str = self._fasta.fetch_string(self.chrom,
                                                 pos, self._GC_model_len-1)
            gc_idx = fasta_str.count('G')+fasta_str.count('C')
            return self._model[gc_idx]
        else:
            return 1
