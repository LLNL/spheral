from math import *
from Spheral import Vector3d
from Spheral import Tensor3d
from Spheral import SymTensor3d

from Spheral import vector_of_int, vector_of_double, vector_of_SymTensor3d, vector_of_vector_of_double

import mpi
rank = mpi.rank
procs = mpi.procs

positions = []

nshell = 1000

if (nshell>procs and procs>1):
    npp = nshell/(procs - 1)
    p = 0
    if(rank>0 and rank*npp!=nshell):
        imax = rank*npp+1
        for i in xrange(1,imax):
            h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
            if(i>1 and i<nshell):
                p = (p+3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h))%(2.0*pi)
    rankmin = rank*npp+1
    rankmax = ((rank+1)*npp + 1) if (rank != procs -1) else (nshell + 1)
    for i in xrange(rankmin, rankmax):
        h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
        t = acos(h)
        if(i>1 and i<nshell):
            p = (p+3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h))%(2.0*pi)
        else:
            p = 0
        x = sin(t)*cos(p)
        y = sin(t)*sin(p)
        z = cos(t)
        positions.append([x,y,z,rank])
else:
    # let rank 0 do all the work
    p = 0
    if (mpi.rank == 0):
        for i in xrange(1, nshell+1):
            h = -1.0+(2.0*(i-1.0)/(nshell-1.0))
            t = acos(h)
            if(i>1 and i<nshell):
                p = (p+3.8/sqrt(nshell)*1.0/sqrt(1.0-h*h))%(2.0*pi)
            else:
                p = 0
            x = sin(t)*cos(p)
            y = sin(t)*sin(p)
            z = cos(t)
            positions.append([x,y,z,0])

positions = mpi.allreduce(positions,mpi.SUM)
#print "x\ty\tz\trho\n"
for i in xrange(len(positions)):
    print "%f\t%f\t%f\t%f"%(positions[i][0],positions[i][1],positions[i][2],positions[i][3])
