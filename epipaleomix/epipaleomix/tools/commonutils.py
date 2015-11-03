import pysam
import re
import os


# def merge_dics(x, y):
#     '''Given two dicts, merge them into a new dict as a shallow copy.'''
#     z = x.copy()
#     z.update(y)
#     return z

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

def _read_mappa(args):
    with open(args.MappabilityPath, 'r') as mappafile:
        for line in mappafile:
            chrom, start, end, score = _unpack_mappa(*(re.split(r'\s+', line.rstrip())))
            yield chrom, start, end, score


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


class GC_correction(object):
    ''' class doc '''
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


# class _GC_correction(object):
#     ''' This is for individual read length '''
#     def __init__(self):
#         if self.arg.GCmodel:
#             half_jump = 4  ## value comes from the 'resolution' used in epaleomix.py # resolution/2
#             self._fasta_dat = Cache(self.arg.FastaPath)
#             with open(self.arg.GCmodel, 'r') as f:
#                 next(f)  # do not need the header
#                 lst = []
#                 self._models_dic = {}
#                 prev_count = ''
#                 for line in f:
#                     curr_count, ratio = self._unpackgc(*re.split(r'\s+', line.rstrip()))

#                     if curr_count != prev_count and lst:
#                         [self._updatedic(key, lst) for key in xrange(prev_count-half_jump,
#                                                                        prev_count+half_jump+1)]
#                         lst=[]
#                     lst.append(ratio)
#                     prev_count = curr_count

#             if lst:
#                 [self._updatedic(key, lst) for key in xrange(prev_count-half_jump,
#                                                              prev_count+half_jump+1)]

#     def _unpackgc(self, curr_count, content, ratio):
#         return int(curr_count), float(ratio)

#     def _updatedic(self, key, lst):
#         self._models_dic[key] = lst

#     def _get_gc_corr_dep(self, record):
#         if self.arg.GCmodel:
#             try:
#                 ## very few reads are aoutside min max length, they are given no enrichment
#                 model = self._models_dic[record.alen]
#             except KeyError:
#                 return 1
#             modellength = len(model)
#             pos = record.aend-1-modellength if record.is_reverse else record.pos
#             # pos = record.aend-1 if record.is_reverse else record.pos
#             fasta_str = self._fasta_dat.fetch_string(self.chrom,
#                                                      pos, modellength-1)
#             gc_idx = fasta_str.count('G')+fasta_str.count('C')
#             return model[gc_idx]

### add this to epaleomix.py if individual readlength to be used

# def concat_gcsubnodes(nodecls, d_bam, gcwindows, subn=()):
#     return [nodecls(d_bam, subnodes=[GccorrectNode(d_bam, rl, subnodes=subn) for rl in gcwindows])]
                
# def calc_gcmodel(d_bam):
#     rlmin, rlmax = getdequelen(d_bam)
#     if d_bam.opts['GCcorrect'].get('Enabled', False):
#         chromused_coerce_to_string(d_bam)
#         checkmappabilitychrom.main([d_bam.prefix.get('--MappabilityPath', MakefileError),
#                                     d_bam.opts['GCcorrect'].get('ChromUsed', MakefileError)])

#         resolution = 9
#         return concat_gcsubnodes(CreateGCModelNode, d_bam,
#                                  xrange(rlmin, rlmax+resolution, resolution))
#     return []
