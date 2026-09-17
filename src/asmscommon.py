import re
import numpy as np
import sys

def parse_clusters(clusterfn):
    """
    returns 2 dictionaries:
    nr[clid] -> number of reads it contains    
    cl[rname] -> clid
    """
    cl = {}
    k = -1
    f = open(clusterfn)
    ## parse header
    while True:
        l = f.readline()
        if ('' == l): break
        l = l.strip()
        m = re.match("^#", l)
        ## expect some info line at the start of the file
        if ( None == m ):
            break
        else:
            m = re.match("#@K:([0-9]+)",l)
            if ( m != None ):
                k = int(m.groups(0)[0])
    if ( -1 == k ):
        print(f"warning: #@K missing from header", file=sys.stderr)
        k=2
        ##sys.exit(1)
    ## parse clusters
    ## nr contains the number of reads in each cluster
    nr = np.full(k, fill_value=0)
    while True:
        if ( '' == l ): break
        ## XXX TODO: only write #@ in the header of the file,
        ## not later
        if (l[0:2] == '#@'):
            l = f.readline()
            print(l)
            continue
        l = l.strip()
        fields = l.split()
        print(fields)
        cln = int(fields[1])
        rname = fields[0]
        l = f.readline()
        cl[rname] = cln
        nr[cln] = nr[cln] + 1
    f.close()
    return (nr, cl)

def gmatrix_of_file(matrixfn):
    """ 
    read matrix into hash tables 
    (which are better for sparse representation)
    """
    if ( '-' == matrixfn ):
        f = sys.stdin
    else:
        f = open(matrixfn, 'r')
    """
    drnames : idx -> rname
    dgpos: idx -> gpos
    dstate: (ridx, gposidx) -> meth
    maxridx: max read index
    maxgposidx: max gpos index
    """
    drnames={}; dgpos={}; dstate={}
    nlines=0; 
    for line in f:
        if ('#' == line[0]) : continue
        line = line.strip()
        fields = line.split()
        """fields: ridx, rname, cidx, gpos, state"""
        ridx = int(fields[0]); cidx = int(fields[2])
        rname = fields[1]; gpos= int(fields[3]); state=int(fields[4])
        drnames[ridx] = rname
        dgpos[cidx] = gpos
        if ( (ridx, cidx) in dstate ):
            print(f"malformed: element ({rname}, {gpos}) appears twice",
                  file = sys.stderr)
            sys.exit(1)
        else:
            dstate[(ridx, cidx)] = state
        nlines+=1
    maxridx = max(list(drnames.keys()))
    maxgposidx = max(list(dgpos.keys()))
    return(drnames, dgpos, dstate, maxridx, maxgposidx)
