from TrampolineGenerator import *

class PyBoundary(TrampolineGenerator):

    def __init__(self):
        TrampolineGenerator.__init__(self)
        self.includes = ["Geometry/GeomPlane.hh",
                         "NodeList/NodeList.hh",
                         "Field/Field.hh",
                         "DataBase/DataBase.hh"]
        self.namespaces = ["Spheral", "BoundarySpace"]
        self.templates = ["Dimension"]
        self.preamble = """
using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
"""
        return

    def setAllGhostNodes(self,
                         args = [("DataBase<Dimension>&", "dataBase")],
                         doc = "Recreate ghost nodes for this Boundary for all NodeLists in the DataBase."):
        return

    def setAllViolationNodes(self,
                             args = [("DataBase<Dimension>&", "dataBase")]):
        return

    def cullGhostNodes(self,
                       args = [("const FieldList<Dimension, int>&", "flagSet"),
                               ("FieldList<Dimension, int>&", "old2newIndexMap"),
                               ("std::vector<int>&", "numNodesRemoved")]):
        return

    def setGhostNodes(self,
                      args = [("NodeList<Dimension>&", "nodeList")],
                      pure = True):
        return

    def updateGhostNodes(self,
                         args = [("NodeList<Dimension>&", "nodeList")],
                         pure = True):
        return

    def applyGhostBoundary1(self,
                            args = [("Field<Dimension, int>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return

    def applyGhostBoundary2(self,
                            args = [("Field<Dimension, Scalar>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return

    def applyGhostBoundary3(self,
                            args = [("Field<Dimension, Vector>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return

    def applyGhostBoundary4(self,
                            args = [("Field<Dimension, Tensor>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return

    def applyGhostBoundary5(self,
                            args = [("Field<Dimension, SymTensor>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return

    def applyGhostBoundary6(self,
                            args = [("Field<Dimension, ThirdRankTensor>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return

    def applyGhostBoundary7(self,
                            args = [("Field<Dimension, std::vector<Scalar>>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return

    def setViolationNodes(self,
                          args = [("NodeList<Dimension>&", "nodeList")],
                          pure = True):
        return

    def updateViolationNodes(self,
                             args = [("NodeList<Dimension>&", "nodeList")],
                             pure = True):
        return

    def enforceBoundary1(self,
                         args = [("Field<Dimension, int>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return

    def enforceBoundary2(self,
                         args = [("Field<Dimension, Scalar>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return

    def enforceBoundary3(self,
                         args = [("Field<Dimension, Vector>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return

    def enforceBoundary4(self,
                         args = [("Field<Dimension, Tensor>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return

    def enforceBoundary5(self,
                         args = [("Field<Dimension, SymTensor>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return

    def enforceBoundary6(self,
                         args = [("Field<Dimension, ThirdRankTensor>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return

    def initializeProblemStartup(self):
        return

    def finalizeGhostBoundary(self,
                              const = True):
        return

    def reset(self,
              args = [("const DataBase<Dimension>&", "dataBase")]):
        return

    def numGhostNodes(self,
                      const = True):
        return "int"

    def clip(self,
             args = [("Vector&", "xmin"),
                     ("Vector&", "xmax")],
             const = True):
        return

# if __name__ == "__main__":
#     generateTrampoline(PyBoundary())
#     generateConcreteTrampoline(PyBoundary())
#     generateBindingFunction(PyBoundary())
