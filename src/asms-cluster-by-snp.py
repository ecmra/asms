#!/usr/bin/env python3

import argparse
import pysam
import sys

"""

read alignment file overlapping the position of the variant
print reads which are ref for that variant and reads which are
alt using the same notations of the asms clusters.

test with:

src/asms-cluster-by-snp.py  /home/ecmra/data/ashkenazi-trio-chr15/son_HG002.chr15.bam chr15:17009256:C:T
"""


def parse_position(s):
    """
    parse chr15:17009156:C:T
    """
    (contig, pos, ref, alt) = s.split(':')
    return (contig, int(pos), ref, alt)
    

parser = argparse.ArgumentParser(description="cluster reads according to snps")
parser.add_argument("bamfn", help="BAM file")
parser.add_argument("var", help="variant description contig:gcoord:ref:alt")
args = parser.parse_args()
bamfn = args.bamfn
var = args.var

with pysam.AlignmentFile(bamfn, "rb") as bam:
    (contig, pos1, ref, alt) = parse_position(var)
    pos = pos1 - 1  # Convert to 0-based
    reads = bam.fetch(contig, pos, pos + 1)
    for read in reads:
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
        if nuc == ref:
            print(f"{rname}\t0\t{1-perror}")
            continue
        if nuc == alt :
            print(f"{rname}\t1\t{1-perror}")
            continue
    
