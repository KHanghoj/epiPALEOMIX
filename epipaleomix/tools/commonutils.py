import pysam
import re
import os


# def merge_dics(x, y):
#     '''Given two dicts, merge them into a new dict as a shallow copy.'''
#     z = x.copy()
#     z.update(y)
#     return z

READ_NUCL = {'non': lambda x: open(x, 'r'), 'gz': lambda x: gzip.open(x, 'rb') }
def _read_bed_nucl(args):
    ''' For EM-based nuclesome calling. takes nucleosome calling as input'''
    assert os.path.exists(args.nuclcalls), "No bedfile available"    
    bed_reader = READ_NUCL['gz'] if re.search(r".gz$", args.nuclcalls) else  READ_NUCL['non']
    if args.nuclcalls:
        with bed_reader(args.nuclcalls) as bedfile:
            line = next(bedfile, None) # need to check if bed file contain a header
            if not re.search(r"^#",line):
                chrom, start, end, _, score, bedc = re.split(r'\s+', line.rstrip())
                if float(score) > 0:
                    yield str(chrom), int(start), int(end)

            for line in bedfile:
                chrom, start, end, _, score, bedc = re.split(r'\s+', line.rstrip())
                if float(score) > 0:
                    yield str(chrom), int(start), int(end)

def unpack_bed_nucl(chrom, start, end, *rest):
    return chrom, start, end
                    
def read_bed_nucl(args):
    ''' For EM-based nuclesome calling. takes nucleosome calling as input'''
    assert os.path.exists(args.nuclcalls), "No bedfile available"    
    bed_reader = READ_NUCL['gz'] if re.search(r".gz$", args.nuclcalls) else  READ_NUCL['non']
    if args.nuclcalls:
        with bed_reader(args.nuclcalls) as bedfile:
            line = next(bedfile, None) # need to check if bed file contain a header
            if not re.search(r"^#",line):
                chrom, start, end = unpack_bed_nucl(*re.split(r'\s+', line.rstrip()))
                yield str(chrom), int(start), int(end)

            for line in bedfile:
                chrom, start, end = unpack_bed_nucl(*re.split(r'\s+', line.rstrip()))
                yield str(chrom), int(start), int(end)


def check_path(temp_dir):
    if not os.path.exists(temp_dir):
        try:
            os.makedirs(temp_dir)
        except OSError, error:
            print_err("ERROR: Could not create temp root:\n\t%s" % (error,))
            return 1


def unpack(chrom, start, end, bedcoord, *rest):
    return chrom, start, end, bedcoord, rest


def unpack_mappa(chrom, start, end, *rest):
    return chrom, start, end, rest


def _unpack_mappa(chrom, start, end, score):
    return chrom, int(start), int(end), float(score)


def read_bed(args):
    if args.bed:
        with open(args.bed, 'r') as bedfile:
            for line in bedfile:
                chrom, start, end, bedc, rest = unpack(*(re.split(r'\s+', line.rstrip())))
                yield str(chrom), int(start), int(end), str(bedc)


def read_mappa(args):
    with open(args.MappabilityPath, 'r') as mappafile:
        for line in mappafile:
            chrom, start, end, rest = unpack_mappa(*(re.split(r'\s+', line.rstrip())))
            yield str(chrom), int(start), int(end), float(rest[-1])

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
            # fetch is 0-based. if 20000000 to 20000001. you get the 20000001th base
            self._fasta_str = self._fasta.fetch(chrom,
                                                start=start, end=self._end)
            self._last_start = start
            self._last_chrom = chrom
        self._actual_pos = start-self._last_start
        return self._fasta_str[self._actual_pos:self._actual_pos+nbases]

    def closefile(self):
        ''' docstring '''
        return self._fasta.close()


class _GC_correction(object):
    ''' This is the original GC_correction '''
    def __init__(self):
        if self.arg.GCmodel:
            self._fasta_dat = Cache(self.arg.FastaPath)
            with open(self.arg.GCmodel, 'r') as f:
                next(f)  # do not need the header
                self._model = [float(line.rstrip('\n').split('\t')[-1])
                               for line in f]
                self._GC_model_len = len(self._model)

    def _get_gc_corr_dep(self, record):
        if self.arg.GCmodel:
            pos = record.aend-self._GC_model_len if record.is_reverse else record.pos
            fasta_str = self._fasta_dat.fetch_string(self.chrom,
                                                     pos, self._GC_model_len-1)
            gc_idx = fasta_str.count('G')+fasta_str.count('C')
            return self._model[gc_idx]
        else:
            return 1


class GC_correction(object):
    ''' This is for individual read length. Just Testing '''
    def __init__(self):
        if self.arg.GCmodel:
            half_jump = 2  ## value comes from the 'resolution' used in epaleomix.py # resolution/2
            self._fasta_dat = Cache(self.arg.FastaPath)
            with open(self.arg.GCmodel, 'r') as f:
                next(f)  # do not need the header
                lst = []
                self._models_dic = {}
                prev_count = ''
                for line in f:
                    curr_count, ratio = self._unpackgc(*re.split(r'\s+', line.rstrip()))

                    if curr_count != prev_count and lst:
                        for key in range(prev_count-half_jump, prev_count+half_jump+1):
                            self._models_dic[key] = lst

                        lst=[]
                    lst.append(ratio)
                    prev_count = curr_count

            if lst:
                for key in range(prev_count-half_jump, prev_count+half_jump+1):
                    self._models_dic[key] = lst

    def _unpackgc(self, curr_count, content, ratio):
        return int(curr_count), float(ratio)

    def _get_gc_corr_dep(self, record):
        if self.arg.GCmodel:
            try:
                model = self._models_dic[record.alen]
            except KeyError:
                if record.alen > max(self._models_dic):
                    model = self._models_dic[max(self._models_dic)]
                else:
                    model = self._models_dic[min(self._models_dic)]
            modellength = len(model)-1 # minus 1 as contains both 0 and max i.e. 50 is 51 long
            pos = record.aend-1-modellength if record.is_reverse else record.pos
            fasta_str = self._fasta_dat.fetch_string(self.chrom,
                                                     pos, modellength)
            gc_idx = fasta_str.count('G')+fasta_str.count('C')
            return model[gc_idx]
        else:
            return 1.0
