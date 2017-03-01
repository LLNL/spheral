import Spheral
import ParMETISDistributeNodes

#-------------------------------------------------------------------------------
# Domain decompose (1d method).
#-------------------------------------------------------------------------------
def distributeNodes1d(listOfNodeTuples):
    return ParMETISDistributeNodes.distributeNodesGeneric(listOfNodeTuples,
                                                          Spheral.DataBase1d,
                                                          Spheral.ScalarField1d,
                                                          Spheral.NestedGridRedistributeNodes1d)

#-------------------------------------------------------------------------------
# Domain decompose (2d method).
#-------------------------------------------------------------------------------
def distributeNodes2d(listOfNodeTuples):
    return ParMETISDistributeNodes.distributeNodesGeneric(listOfNodeTuples,
                                                          Spheral.DataBase2d,
                                                          Spheral.ScalarField2d,
                                                          Spheral.NestedGridRedistributeNodes2d)

#-------------------------------------------------------------------------------
# Domain decompose (3d method).
#-------------------------------------------------------------------------------
def distributeNodes3d(listOfNodeTuples):
    return ParMETISDistributeNodes.distributeNodesGeneric(listOfNodeTuples,
                                                          Spheral.DataBase3d,
                                                          Spheral.ScalarField3d,
                                                          Spheral.NestedGridRedistributeNodes3d)
