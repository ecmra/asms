#!/usr/bin/env python3

import argparse
import numpy as np
import scipy
import sys



"""
read two clusters
for all the reads which coincide across the clusters
fill an array with the counts of (0,0),(0,1),(1,0),(1,1)
and compute the pvalue of a chi-squared test on it"
"""

parser = argparse.ArgumentParser(description="")
parser.add_argument("listfn")


"""
test with:
src/asms-cmp-clusters.py ../asms-tests/out/cluster-by-snp-test/matched-clusters.txt 

matched-clusters is a 2 columns file.
first column: snp file
second column: meth cluster file

"""

def read_cluster(fn):
    d = {}
    with open(fn, 'r') as f:
        for line  in f:
            if line[0] == '#':
                continue
            fields = line.strip().split('\t')
            d[fields[0]] = int(fields[1])
    return d

def make_counts_array(dsnp, dmeth):
    ck = set(dsnp.keys()).intersection(set(dmeth.keys()))
    a = np.zeros((2,2), dtype=int)
    for k in ck:
        a[dsnp[k],dmeth[k]]+=1
    return a


def chisquare_test(a):
    a = a.reshape(4)
    t = scipy.stats.chisquare(a, ddof=1)
    return (a, t.statistic, t.pvalue)


if __name__ == '__main__':
    args = parser.parse_args()
    listfn = args.listfn
    listf = open(listfn ,'r')
    for line in listf:
        fields = line.strip().split('\t')
        fnsnp = fields[0]; fnmeth = fields[1]
        dsnp = read_cluster(fnsnp)
        dmeth = read_cluster(fnmeth)
        a = make_counts_array(dsnp,dmeth)
        pop0snp = a[0,0] + a[0,1]; pop1snp=a[1,0] + a[1,1]
        pop0meth = a[0,0]+a[1,0]; pop1meth=a[0,1]+a[1,1]
        if pop0snp < pop1snp:
            af = pop0snp/(pop0snp + pop1snp)
        else:
            af = pop1snp/(pop0snp + pop1snp)
        if pop0meth < pop1meth:
            mratio = pop0meth/(pop0meth + pop1meth)
        else:
            mratio = pop1meth/(pop0meth + pop1meth)
        (a,stat, pval) = chisquare_test(a)
        print(f"{fnsnp}\t{fnmeth}\t{a}\t{stat:.3g}\t{pval:.3g}\t{af:.3f}\t{mratio:.3f}")


