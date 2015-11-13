#!/usr/bin/python
#
# Copyright (c) 2012 Mikkel Schubert <MSchubert@snm.ku.dk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
import os
import sys
import argparse
import subprocess

from pypeline.common.fileutils import missing_executables


def _collect_positions(handle):
    positions = set()
    for line in handle:
        if not line.startswith("#"):
            fields = line.split("\t", 6)

            if "," in fields[4]:
                positions.add((fields[0], int(fields[1])))
    positions = list(positions)
    positions.sort()

    return positions


def print_status(index, contig, position, npositions, end="\r"):
    tmpl = " - Collected pileup %i of %i (%02.2f%% done): %s:%s ...      %s"
    sys.stderr.write(tmpl % (index,
                             npositions,
                             (100.0 * index) / npositions,
                             contig,
                             position,
                             end))


def main(argv):
    parser = argparse.ArgumentParser(prog="paleomix create_pileup")
    parser.add_argument("output", help="BGZipped pileup file.")
    parser.add_argument("mpileup_args", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)

    missing = missing_executables(("samtools", "tabix", "bgzip"))
    if missing:
        sys.stderr.write("ERROR: Required executables are missing:\n")
        sys.stderr.write("    - %s\n" % "\n\t- ".join(missing))
        return 1

    with open(args.output, "w") as handle:
        bgzip = subprocess.Popen("bgzip",
                                 stdin  = subprocess.PIPE,
                                 stdout = handle)

        # While samtools mpileup has an option for specifying a list of
        # positions (-l), this requires traversing the entire file, and may
        # not calculate the BAQ. Given the low number of expected sites,
        # individual calls for each position are significantly faster.
        sys.stderr.write("Reading VCF from STDIN ...\n")
        positions = _collect_positions(sys.stdin)
        npositions = len(positions)
        sys.stderr.write(" - Read %i candidate positions ...\n"
                         % (npositions,))

        positions_file = args.output + ".positions"
        with open(positions_file, "w") as handle:
            for (contig, position) in positions:
                handle.write("%s\t%s\n" % (contig, position))
        sys.stderr.write(" - Wrote positions to '%s' ...\n"
                         % (positions_file,))

        sys.stderr.write("Collecting pileups:\n")
        call = ["samtools", "mpileup", "-R", "-l", positions_file]
        call.extend(args.mpileup_args)
        proc = subprocess.Popen(call,
                                stdout=subprocess.PIPE,
                                close_fds=True)

        line = "NA\tNA\t"
        index = npositions = 1
        for (index, line) in enumerate(proc.stdout, start=1):
            if (index - 1) % 100 == 0:
                contig, position, _ = line.split("\t", 2)
                print_status(index, contig, position, npositions)
            bgzip.stdin.write(line)

        contig, position, _ = line.split("\t", 2)
        print_status(index, contig, position, npositions, end="\n")

        if proc.wait():
            sys.stderr.write("ERROR: Error running samtools, return-code=%i:\n"
                             % proc.wait())
            sys.stderr.write("\t- Command: %s\n" % " ".join(call))
            return 1

        bgzip.stdin.close()
        if bgzip.wait():
            sys.stderr.write("ERROR: Error running bgzip, return-code %i\n"
                             % bgzip.wait())

    sys.stderr.write(" - Cleaning up ...")
    os.remove(positions_file)

    subprocess.check_call(["tabix", "-b", "2", "-e", "2", args.output])
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
