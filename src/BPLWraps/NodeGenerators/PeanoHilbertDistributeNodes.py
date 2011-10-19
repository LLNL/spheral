import Spheral
import distributeNodesGeneric

#-------------------------------------------------------------------------------
# Domain decompose using PeanoHilbert ordering (1d method).
#-------------------------------------------------------------------------------
def distributeNodes1d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase1d,
                                                  Spheral.globalNodeIDs,
                                                  Spheral.PeanoHilbertOrderRedistributeNodes1d)

#-------------------------------------------------------------------------------
# Domain decompose using PeanoHilbert ordering (2d method).
#-------------------------------------------------------------------------------
def distributeNodes2d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase2d,
                                                  Spheral.globalNodeIDs,
                                                  Spheral.PeanoHilbertOrderRedistributeNodes2d)

#-------------------------------------------------------------------------------
# Domain decompose using PeanoHilbert ordering (3d method).
#-------------------------------------------------------------------------------
def distributeNodes3d(*listOfNodeTuples):
    distributeNodesGeneric.distributeNodesGeneric(listOfNodeTuples,
                                                  Spheral.DataBase3d,
                                                  Spheral.globalNodeIDs,
                                                  Spheral.PeanoHilbertOrderRedistributeNodes3d)

