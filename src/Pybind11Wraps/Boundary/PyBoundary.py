import sys
sys.path.append("..")
from TrampolineGenerator import TrampolineGenerator

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

    def setGhostNodes(self,
                      args = [("NodeList<Dimension>&", "nodeList")]):
        return "void"

if __name__ == "__main__":
    PyBoundary()()
