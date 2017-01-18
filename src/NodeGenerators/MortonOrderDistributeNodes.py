import Spheral
import distributeNodesGeneric

#-------------------------------------------------------------------------------
# Domain decompose using Morton ordering (1d method).
#-------------------------------------------------------------------------------
def distributeNodes1d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase1d,
                                                  Spheral.globalNodeIDsAll1d,
                                                  Spheral.MortonOrderRedistributeNodes1d)

#-------------------------------------------------------------------------------
# Domain decompose using Morton ordering (2d method).
#-------------------------------------------------------------------------------
def distributeNodes2d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase2d,
                                                  Spheral.globalNodeIDsAll2d,
                                                  Spheral.MortonOrderRedistributeNodes2d)

#-------------------------------------------------------------------------------
# Domain decompose using Morton ordering (3d method).
#-------------------------------------------------------------------------------
def distributeNodes3d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase3d,
                                                  Spheral.globalNodeIDsAll3d,
                                                  Spheral.MortonOrderRedistributeNodes3d)

