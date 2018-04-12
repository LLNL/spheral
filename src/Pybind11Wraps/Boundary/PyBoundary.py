from TrampolineGenerator import *

class PyBoundary(TrampolineGenerator):

    def __init__(self):
        TrampolineGenerator.__init__(self,
                                     includes = ["Geometry/GeomPlane.hh",
                                                 "NodeList/NodeList.hh",
                                                 "Field/Field.hh",
                                                 "DataBase/DataBase.hh"],
                                     namespaces = ["Spheral", "BoundarySpace"],
                                     templates = ["Dimension"],
                                     preamble = """
using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
""")
        return

    def setAllGhostNodes(self,
                         args = [("DataBase<Dimension>&", "dataBase")],
                         doc = "Recreate ghost nodes for this Boundary for all NodeLists in the DataBase."):
        return "void"

    def setAllViolationNodes(self,
                             args = [("DataBase<Dimension>&", "dataBase")]):
        return "void"

    def cullGhostNodes(self,
                       args = [("const FieldList<Dimension, int>&", "flagSet"),
                               ("FieldList<Dimension, int>&", "old2newIndexMap"),
                               ("std::vector<int>&", "numNodesRemoved")]):
        return "void"

    def setGhostNodes(self,
                      args = [("NodeList<Dimension>&", "nodeList")],
                      pure = True):
        return "void"

    def updateGhostNodes(self,
                         args = [("NodeList<Dimension>&", "nodeList")],
                         pure = True):
        return "void"

    def applyGhostBoundary1(self,
                            args = [("Field<Dimension, int>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return "void"

    def applyGhostBoundary2(self,
                            args = [("Field<Dimension, Scalar>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return "void"

    def applyGhostBoundary3(self,
                            args = [("Field<Dimension, Vector>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return "void"

    def applyGhostBoundary4(self,
                            args = [("Field<Dimension, Tensor>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return "void"

    def applyGhostBoundary5(self,
                            args = [("Field<Dimension, SymTensor>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return "void"

    def applyGhostBoundary6(self,
                            args = [("Field<Dimension, ThirdRankTensor>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return "void"

    def applyGhostBoundary7(self,
                            args = [("Field<Dimension, std::vector<Scalar>>&", "field")],
                            const = True,
                            pure = True,
                            name = "applyGhostBoundary"):
        return "void"

    def setViolationNodes(self,
                          args = [("NodeList<Dimension>&", "nodeList")],
                          pure = True):
        return "void"

    def updateViolationNodes(self,
                             args = [("NodeList<Dimension>&", "nodeList")],
                             pure = True):
        return "void"

    def enforceBoundary1(self,
                         args = [("Field<Dimension, int>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary2(self,
                         args = [("Field<Dimension, Scalar>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary3(self,
                         args = [("Field<Dimension, Vector>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary4(self,
                         args = [("Field<Dimension, Tensor>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary5(self,
                         args = [("Field<Dimension, SymTensor>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary6(self,
                         args = [("Field<Dimension, ThirdRankTensor>&", "field")],
                         const = True,
                         pure = True,
                         name = "enforceBoundary"):
        return "void"

    def initializeProblemStartup(self):
        return "void"

    def finalizeGhostBoundary(self,
                              const = True):
        return "void"

    def reset(self,
              args = [("const DataBase<Dimension>&", "dataBase")]):
        return "void"

    def numGhostNodes(self,
                      const = True):
        return "int"

    def clip(self,
             args = [("Vector&", "xmin"),
                     ("Vector&", "xmax")],
             const = True):
        return "void"

# if __name__ == "__main__":
#     generateTrampoline(PyBoundary())
#     generateConcreteTrampoline(PyBoundary())
#     generateBindingFunction(PyBoundary())
