import Spheral
import distributeNodesGeneric

#-------------------------------------------------------------------------------
# Domain decompose using PeanoHilbert ordering (1d method).
#-------------------------------------------------------------------------------
def distributeNodes1d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase1d,
                                                  Spheral.globalNodeIDsAll1d,
                                                  Spheral.PeanoHilbertOrderRedistributeNodes1d)

#-------------------------------------------------------------------------------
# Domain decompose using PeanoHilbert ordering (2d method).
#-------------------------------------------------------------------------------
def distributeNodes2d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase2d,
                                                  Spheral.globalNodeIDsAll2d,
                                                  Spheral.PeanoHilbertOrderRedistributeNodes2d)

#-------------------------------------------------------------------------------
# Domain decompose using PeanoHilbert ordering (3d method).
#-------------------------------------------------------------------------------
def distributeNodes3d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase3d,
                                                  Spheral.globalNodeIDsAll3d,
                                                  Spheral.PeanoHilbertOrderRedistributeNodes3d)

