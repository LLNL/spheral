from TrampolineGenerator import *

class PyPlanarBoundary(TrampolineGenerator):

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

    def setGhostNodes(self,
                      args = [("NodeList<Dimension>&", "nodeList")]):
        return "void"

    def updateGhostNodes(self,
                           args = [("NodeList<Dimension>&", "nodeList")]):
        return "void"

    def setViolationNodes(self,
                          args = [("NodeList<Dimension>&", "nodeList")]):
        return "void"

    def updateViolationNodes(self,
                             args = [("NodeList<Dimension>&", "nodeList")]):
        return "void"

    def valid(self,
              const = True):
        return "bool"
