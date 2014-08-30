import Spheral
import distributeNodesGeneric

# Import mpi.
import mpi
rank = mpi.rank
procs = mpi.procs

#-------------------------------------------------------------------------------
# Domain decompose using Voronoi ordering (1d method).
#-------------------------------------------------------------------------------
def distributeNodes1d(*listOfNodeTuples):
    if procs > 1:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase1d,
                                                      Spheral.globalNodeIDsAll1d,
                                                      Spheral.VoronoiRedistributeNodes1d)
    else:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase1d,
                                                      Spheral.globalNodeIDsAll1d,
                                                      None)

#-------------------------------------------------------------------------------
# Domain decompose using Voronoi ordering (2d method).
#-------------------------------------------------------------------------------
def distributeNodes2d(*listOfNodeTuples):
    if procs > 1:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase2d,
                                                      Spheral.globalNodeIDsAll2d,
                                                      Spheral.VoronoiRedistributeNodes2d)
    else:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase2d,
                                                      Spheral.globalNodeIDsAll2d,
                                                      None)

#-------------------------------------------------------------------------------
# Domain decompose using Voronoi ordering (3d method).
#-------------------------------------------------------------------------------
def distributeNodes3d(*listOfNodeTuples):
    if procs > 1:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase3d,
                                                      Spheral.globalNodeIDsAll3d,
                                                      Spheral.VoronoiRedistributeNodes3d)
    else:
        distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                      Spheral.DataBase3d,
                                                      Spheral.globalNodeIDsAll3d,
                                                      None)
