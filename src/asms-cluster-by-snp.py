#!/usr/bin/env python3

import argparse
import pysam
import sys
import time

"""
read alignment file overlapping the position of the variant
print reads which are ref for that variant and reads which are
alt using the same notations of the asms clusters.

test with:
src/asms-cluster-by-snp.py ../asms-tests/out/cluster-by-snp-test/bam-snps.txt  ../asms-tests/out/cluster-by-snp-test/test
"""

def parse_position(s):
    """
    parse chr15:17009156:C:T
    """
    (contig, pos, ref, alt) = s.split(':')
    return (contig, int(pos), ref, alt)

parser = argparse.ArgumentParser(description="cluster reads according to snps")
parser.add_argument("bamvarfn", help="BAM files and SNPs (contig:gcoord:ref:alt)")
parser.add_argument("outprefix", help="prefix for output files")
args = parser.parse_args()
bamvarfn = args.bamvarfn
outprefix = args.outprefix

def cluster_by_snp(bamfn, var, outf):
    with pysam.AlignmentFile(bamfn, "rb") as bam:
        (contig, pos1, ref, alt) = parse_position(var)
        print((contig, pos1, ref, alt), file=sys.stderr)
        pos0 = pos1 - 1  # Convert to 0-based
        reads = bam.fetch(contig, pos0, pos0 + 1)
        for read in reads:
            rname = read.query_name
            flag = read.flag
            """
            BAM_FUNMAP  4
            BAM_FSECONDARY 256
            BAM_FSUPPLEMENTARY 2048
            BAM_FQCFAIL 512
            BAM_FDUP    1024
            """
            qual = read.mapping_quality
            if flag & (4 | 256 | 2048 | 512 | 1024):
                continue
            if qual < 10:
                continue
            rpos_gpos = read.get_aligned_pairs(matches_only=True)
            try:
                idx = [e[1] for e in rpos_gpos].index(pos0)
            except ValueError:
                continue
            try:
                # return None if there's no associated sequence
                ##nuc = read.get_forward_sequence()[idx]
                ##q   = read.get_forward_qualities()[idx]
                nuc = read.query_sequence[rpos_gpos[idx][0]]
                q =   read.query_qualities[rpos_gpos[idx][0]]
            except TypeError:
                print(f"warning:{rname}:{flag} has no associated sequence",\
                      file=sys.stderr)
                continue
            perror = 10 ** (-q / 10)
            if nuc == ref:
                print(f"{rname}\t0\t{1-perror}\t{nuc}", file = outf)
            elif nuc == alt :
                print(f"{rname}\t1\t{1-perror}\t{nuc}", file = outf)
            

if bamvarfn == "-":
    bamvarf = sys.stdin
else:
    bamvarf = open(bamvarfn, 'r') 


start = time.time()
lc = 0
for line in bamvarf:
    fields = line.strip().split('\t')
    bamfn = fields[0]; var = fields[1]
    svar = var.replace(':','_')
    outfn=f"{outprefix}-{svar}.txt"
    outf = open(outfn, 'w')
    cluster_by_snp(fields[0], fields[1], outf)
    outf.close()
    lc+=1
    if lc%100 == 0:
        elapsed=time.time()-start
        print(f"{elapsed:.3f}:{lc} variants done", file = sys.stderr)

if bamvarfn != "-":
    bamvarf.close()

statsfn = f"{outprefix}-summary.txt"
statsf= open(statsfn, 'w')
print(f"clustered {lc} variant(s)", file=statsf)
elapsed = time.time() - start
print(f"elapsed:{elapsed:.3f}s", file=statsf)
