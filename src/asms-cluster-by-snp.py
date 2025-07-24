#!/usr/bin/env python3

import argparse
import pysam
import sys

"""
test on
/home/ecmra/data/ashkenazi-trio-chr15/son_HG002.chr15.bam
../asms-tests/out/simple-asms-test/son_HG002.chr15-meth-at-snp.bed
chr15	17009156	17009356	0.4727	3	17009256	C>T
"""


def parse_positions(posfn):
    """
    return a list of positions, a list of ref nucleotides
    and a list of alt nucleotides.
    """
    lpos = []
    with open(posfn) as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom = fields[0]
            pos = int(fields[5])
            ref = fields[6].split('>')[0]
            alt = fields[6].split('>')[1]
            lpos.append((chrom, pos, ref, alt))
    return lpos


parser = argparse.ArgumentParser(description="cluster reads according to snps")
parser.add_argument("bamfn", help="BAM file")
parser.add_argument("posfn", help="bed file with positions")
args = parser.parse_args()
bamfn = args.bamfn
posfn = args.posfn

l = parse_positions(posfn)
with pysam.AlignmentFile(bamfn, "rb") as bam:
    sys.stdout.write("rname\tcluster\prob")
    for chrom, pos1 , ref, alt in l:
        pos = pos1 - 1  # Convert to 0-based
        for read in bam.fetch(chrom, pos, pos + 1):
            rname = read.query_name
            flag = read.flag
            ##print(rname)
            rpos_gpos = read.get_aligned_pairs(matches_only=True)
            try:
                idx = [e[1] for e in rpos_gpos].index(pos)
            except ValueError:
                print(f"warning pos {pos1} not in {rname}",\
                      file=sys.stderr)
                continue
            try:
                nuc = read.get_forward_sequence()[idx]
                q = read.get_forward_qualities()[idx]
            except TypeError:
                print(f"warning:{rname}",\
                      file=sys.stderr)
                continue
            perror = 10 ** (-q / 10)
            print(f"{rname}\t{pos1}\t{flag}\t{ref}\t{alt}\t{nuc}\t{perror}")
