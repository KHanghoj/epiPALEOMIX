#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (c) 2014 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
import sys
import random
import argparse
import itertools

from pypeline.common.sampling import \
    reservoir_sampling, \
    weighted_sampling


###############################################################################
###############################################################################

def _build_rev_compl_table():
    table = ["N"] * 256
    for nt_a, nt_b in zip("AC", "TG"):
        table[ord(nt_a)] = nt_b
        table[ord(nt_b)] = nt_a
    return "".join(table)
_COMPL_TABLE = _build_rev_compl_table()


def _reverse_complement(sequence):
    return sequence.translate(_COMPL_TABLE)[::-1]


###############################################################################
###############################################################################

class ErrorModel(object):
    def __init__(self, offset=33):
        self._offset = offset
        self._models = self._construct_error_models()

    def __getitem__(self, quality):
        idx = ord(quality) - self._offset
        return self._models[idx]

    @classmethod
    def _construct_error_models(cls):
        models = []
        for phred in xrange(0, 100):
            errors = {}
            error_rate = min(0.75, 10 ** (phred / -10.0))
            for index, nuc in enumerate("ACGT"):
                choices = "ACGT"
                weights = [error_rate / 3] * 4
                weights[index] = 1.0 - error_rate
                errors[nuc] = weighted_sampling(choices, weights)
            models.append(errors)
        return models


class RandomFASTA(object):
    def __init__(self, filename, indel_rate=None):
        self._errors = ErrorModel()
        self._seqs = self._read_fasta(filename)
        self._indel_rate = indel_rate

        self._n_errors = 0
        self._n_total = 0
        self._n_indels = 0

    def generate_read(self, read_id, insert_size, qual_1, qual_2, pcr1, pcr2):
        insert = self._fetch(insert_size)
        insert_fw = insert + pcr1
        insert_rw = _reverse_complement(insert) + pcr2

        read_fw = self._read_sequence(insert_fw, qual_1)
        read_rw = self._read_sequence(insert_rw, qual_2)

        return {"fw":  (read_fw, qual_1),
                "rw":  (read_rw, qual_2),
                "ins_len": insert_size,
                "id":  read_id}

    def summary(self):
        print "#Errors #TotalBP #ErrorRate #IndelRate"
        print self._n_errors, self._n_total,
        print float(self._n_errors) / self._n_total,
        print float(self._n_indels) / self._n_total

    def _fetch(self, insert_size):
        candidates = [seq for seq in self._seqs if len(seq) >= insert_size]
        if not candidates:
            raise RuntimeError("Inserts longer than all contigs requested")

        sequence = "N"
        while "N" in sequence:
            candidate = random.choice(candidates)
            start = random.randint(0, len(candidate) - insert_size)
            sequence = candidate[start:start + insert_size]

        assert len(sequence) == insert_size, (len(sequence), insert_size)
        return sequence

    def _read_sequence(self, sequence, qualities, padding="A"):
        read = []
        seq_iter = iter(itertools.chain(sequence, itertools.cycle(padding)))
        qual_iter = iter(qualities)
        while len(read) < len(qualities):
            self._n_total += 1

            if self._indel_rate and self._indel_rate >= random.random():
                self._n_indels += 1

                if random.random() < 0.5:
                    # Deletion
                    seq_iter.next()
                    continue
                else:
                    # Insertion
                    observed = random.choice("ACGT")
                    qual_iter.next()
            else:
                base = seq_iter.next()
                quality = qual_iter.next()
                err_dist = self._errors[quality]
                observed = err_dist[base].next()
                if base != observed:
                    self._n_errors += 1

            read.append(observed)

        return "".join(read)

    @classmethod
    def _read_fasta(cls, filename):
        sequences = []
        with open(filename) as handle:
            current = []
            for line in handle:
                if line.startswith(">"):
                    sequences.append("".join(current))
                    print "Reading", line[1:].split()[0]
                    current = []
                else:
                    current.append(line.rstrip().upper())

            if current:
                sequences.append("".join(current))
            print "Done ..."
        return sequences



###############################################################################
###############################################################################

def _read_qualities(mate_1, mate_2):
    def _read_quality(filename):
        with open(filename) as handle:
            for qualities in itertools.islice(handle, 3, None, 4):
                yield qualities.rstrip()

    quals_1 = _read_quality(mate_1)
    quals_2 = _read_quality(mate_2)
    for (qual_1, qual_2) in itertools.izip(quals_1, quals_2):
        yield qual_1, qual_2


def _write_read(args, read, mate, handle):
    tmpl = "@{prefix}_{id}/{mate} read_len={read_len};ins_len={ins_len}\n"
    ins_len = read["ins_len"]
    sequence, qualities = read[mate]
    handle.write(tmpl.format(read_len=len(sequence),
                             prefix=args.read_name_prefix,
                             mate=1 if mate == "fw" else 2,
                             **read))
    handle.write(sequence)
    handle.write("\n+\n")
    handle.write(qualities)
    handle.write("\n")


def parse_args(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('fasta_ref')
    parser.add_argument('fastq_1')
    parser.add_argument('fastq_2')

    parser.add_argument('output_prefix')
    parser.add_argument('--read-name-prefix', default="Read",
                        help="String prefixed to read names [%(default)s]")
    parser.add_argument('--insert-size-mean', default=[], action="append",
                        help="Insert size mean [100]", type=int)
    parser.add_argument('--insert-size-sd',   default=[], action="append",
                        help="Insert size sd [20]", type=int)
    parser.add_argument('--indel-rate', default=None, type=float,
                        help="Rate of indels [default: off]")
    parser.add_argument('--downsample', metavar="N", type=int, default=None,
                        help="Select only N random templates [default: off]")
    parser.add_argument('--pcr1', default=[], action="append")
    parser.add_argument('--pcr2', default=[], action="append")

    args = parser.parse_args(argv)
    assert len(args.pcr1) == len(args.pcr2)
    if not args.pcr1:
        args.pcr1.append("AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG")
        args.pcr2.append("AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT")

    if not args.insert_size_mean:
        args.insert_size_mean.append(100)
    if not args.insert_size_sd:
        args.insert_size_sd.append(20)

    assert len(args.insert_size_mean) == len(args.insert_size_sd)
    args.insert_sizes = zip(args.insert_size_mean, args.insert_size_sd)

    # Avoid accidentially using these ...
    del args.insert_size_mean
    del args.insert_size_sd

    print "Generating reads with"
    for (mean, sd) in args.insert_sizes:
        print "    Mean = %i; SD = %i for ~1/%i of reads" \
            % (mean, sd, len(args.insert_sizes))

    for adapters in (args.pcr1, args.pcr2):
        for idx, adapter in enumerate(adapters):
            final = []
            for nucl in adapter.upper():
                if nucl == "N":
                    nucl = random.choice("ACGT")
                final.append(nucl)
            adapters[idx] = "".join(final)

    for (pcr1, pcr2) in zip(args.pcr1, args.pcr2):
        print "Generating reads with adapters (~1/%i of reads):" \
            % (len(args.pcr1),)
        print "    --pcr1: %s" % pcr1
        print "    --pcr2: %s" % pcr2

    return args


def read_quality_pairs(args):
    # Generate a synthetic read for each sequence in the templates
    read_id = 0
    for qual_1, qual_2 in _read_qualities(args.fastq_1, args.fastq_2):
        yield read_id, qual_1, qual_2
        read_id += 1


def main(argv):
    args = parse_args(argv)
    fasta = RandomFASTA(args.fasta_ref, indel_rate=args.indel_rate)
    adapters = zip(args.pcr1, args.pcr2)

    quality_pairs = read_quality_pairs(args)
    if args.downsample:
        quality_pairs = reservoir_sampling(quality_pairs, args.downsample)

    with open(args.output_prefix + "_1.fastq", "w") as handle_1:
        with open(args.output_prefix + "_2.fastq", "w") as handle_2:
            for read_id, qual_1, qual_2 in quality_pairs:
                mu, sigma = random.choice(args.insert_sizes)
                insert_size = max(len(qual_1), len(qual_2))

                pcr1, pcr2 = random.choice(adapters)

                read = fasta.generate_read(read_id, insert_size,
                                           qual_1, qual_2,
                                           pcr1, pcr2)

                _write_read(args, read, "fw", handle_1)
                _write_read(args, read, "rw", handle_2)

    with open(args.output_prefix + ".adapters", "w") as handle:
        for pcr1, pcr2 in zip(args.pcr1, args.pcr2):
            handle.write("%s\t%s\n" % (pcr1, pcr2))

    fasta.summary()

    return 0


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
