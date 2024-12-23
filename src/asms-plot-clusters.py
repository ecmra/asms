#!/usr/bin/env python3

import argparse
import asmscommon
import datetime
import matplotlib.pylab as plt
from matplotlib.colors import ListedColormap
import numpy as np
import os
import sys




parser = argparse.ArgumentParser(description="plots 2 clusters generated by cvlr")
parser.add_argument("--save", help="save figure")
parser.add_argument("clusterfn",  help="cluster file")
parser.add_argument("matrixfn",  help="matrix file")

args = parser.parse_args()

clusterfn = args.clusterfn
matrixfn = args.matrixfn

if args.save:
    outfile = args.save
    print( f"going to save in {outfile}", file = sys.stderr )

"""
array nr[clid] -> number of reads it contains    
dict cl[rname] -> clid
"""

nr, cl = asmscommon.parse_clusters(clusterfn)


"""
    drnames : idx -> rname
    dgpos: idx -> gpos
    dstate: (ridx, gposidx) -> meth
    maxridx: max read index
    maxgposidx: max gpos index
"""

drnames, dgpos, dstate, maxridx, maxgposidx = asmscommon.gmatrix_of_file(matrixfn)




m0 = np.full((nr[0], maxgposidx + 1), fill_value=-1, dtype = int)
m1 = np.full((nr[1], maxgposidx + 1), fill_value=-1, dtype = int)

idx0=0; idx1=0
for i in range(maxridx+1):
    rname = drnames[i]
    if ( 0 == cl[rname] ):
        v = m0[idx0,:]
        idx0=idx0+1
    elif ( 1 == cl[rname] ):
        v= m1[idx1,:]
        idx1=idx1+1
    else:
        print(f"no cluster for read :{rname}", file=sys.stderr)
        sys.exit(1)
    for j in range(maxgposidx+1):
        if (i,j) in dstate:
            v[j] = dstate[(i,j)]
            
            
cmap = ListedColormap(["gray", "orange", "blue"])

fig, axs = plt.subplots(2, 1, sharex=True)
cmap = ListedColormap(["gray", "orange", "blue"])
axs[0].imshow(m0, cmap=cmap, vmin=-1, vmax=1, aspect='auto')
axs[0].axis('off')
p1 = axs[1].imshow(m1, cmap=cmap, vmin=-1, vmax=1, aspect='auto')
axs[1].axis('off')

cbar = fig.colorbar(p1,  ticks=[-1, 0, 1], orientation='horizontal', aspect=30, extend='neither', spacing='uniform',
                    drawedges=False, boundaries=[-2,-1,0,1], values=[-1,0,1])
labels = ["NA", "0", "1"]
cbar.ax.set_xticks([-1.5, -0.5, 0.5], labels) 

fig.suptitle(f"cluster0:{m0.shape} cluster1:{m1.shape}")

if args.save:
    plt.savefig(args.save)
else:
    plt.show()
