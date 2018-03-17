import sys
sys.path.append("..")
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
                         args = [("DataBase<Dimension>&", "dataBase")]):
        return "void"

    def setAllViolationNodes(self,
                             args = [("DataBase<Dimension>&", "dataBase")]):
        return "void"

    def cullGhostNodes(self,
                       args = [("const FieldList<Dimension, int>", "flagSet"),
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

    def applyGhostBoundary(self,
                           args = [("Field<Dimension, int>&", "field")],
                           const = True,
                           pure = True):
        return "void"

    def applyGhostBoundary(self,
                           args = [("Field<Dimension, Scalar>&", "field")],
                           const = True,
                           pure = True):
        return "void"

    def applyGhostBoundary(self,
                           args = [("Field<Dimension, Vector>&", "field")],
                           const = True,
                           pure = True):
        return "void"

    def applyGhostBoundary(self,
                           args = [("Field<Dimension, Tensor>&", "field")],
                           const = True,
                           pure = True):
        return "void"

    def applyGhostBoundary(self,
                           args = [("Field<Dimension, SymTensor>&", "field")],
                           const = True,
                           pure = True):
        return "void"

    def applyGhostBoundary(self,
                           args = [("Field<Dimension, ThirdRankTensor>&", "field")],
                           const = True,
                           pure = True):
        return "void"

    def applyGhostBoundary(self,
                             args = [("Field<Dimension, std::vector<Scalar>>&", "field")],
                             const = True,
                             pure = True):
        return "void"

    def setViolationNodes(self,
                          args = [("NodeList<Dimension>&", "nodeList")],
                          pure = True):
        return "void"

    def updateViolationNodes(self,
                             args = [("NodeList<Dimension>&", "nodeList")],
                             pure = True):
        return "void"

    def enforceBoundary(self,
                          args = [("Field<Dimension, int>&", "field")],
                          const = True,
                          pure = True):
        return "void"

    def enforceBoundary(self,
                        args = [("Field<Dimension, Scalar>&", "field")],
                        const = True,
                        pure = True):
        return "void"

    def enforceBoundary(self,
                        args = [("Field<Dimension, Vector>&", "field")],
                        const = True,
                        pure = True):
        return "void"

    def enforceBoundary(self,
                        args = [("Field<Dimension, Tensor>&", "field")],
                        const = True,
                        pure = True):
        return "void"

    def enforceBoundary(self,
                        args = [("Field<Dimension, SymTensor>&", "field")],
                        const = True,
                        pure = True):
        return "void"

    def enforceBoundary(self,
                        args = [("Field<Dimension, ThirdRankTensor>&", "field")],
                        const = True,
                        pure = True):
        return "void"

    def initializeProblemStartup(self):
        return "void"

    def finalizeGhostBoundary(self):
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

if __name__ == "__main__":
    generateAbstractTrampoline(PyBoundary())
    generateConcreteTrampoline(PyBoundary())
    generateBindingFunction(PyBoundary())
