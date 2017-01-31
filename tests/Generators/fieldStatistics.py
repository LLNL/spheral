# Return the min, max, avg, and std deviation for a field.

from math import *
import mpi

def fieldStatistics(f):
    f2vals = [x*x for x in f.internalValues()]
    n = mpi.allreduce(f.numInternalElements, mpi.SUM)
    f2avg = mpi.allreduce(sum(f2vals + [0.0]), mpi.SUM)/max(1, n)
    favg = f.sumElements()/max(1, n)
    return f.min(), f.max(), favg, sqrt(abs(f2avg - favg*favg))
