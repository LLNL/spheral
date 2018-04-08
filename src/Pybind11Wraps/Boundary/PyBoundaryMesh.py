from TrampolineGenerator import *

class PyBoundaryMesh(TrampolineGenerator):

    def __init__(self):
        TrampolineGenerator.__init__(self,
                                     includes = ["Geometry/GeomPlane.hh",
                                                 "NodeList/NodeList.hh",
                                                 "Field/Field.hh",
                                                 "DataBase/DataBase.hh",
                                                 "Field/FieldList.hh",
                                                 "Mesh/Mesh.hh"],
                                     namespaces = ["Spheral", "BoundarySpace"],
                                     templates = ["Dimension"],
                                     preamble = """
using Spheral::NodeSpace::NodeList;
using Spheral::FieldSpace::Field;
using Spheral::FieldSpace::FieldList;
using Spheral::DataBaseSpace::DataBase;
using Spheral::MeshSpace::Mesh;
""")
        return

    def enforceBoundary1(self,
                         args = [("std::vector<int>&", "faceField"),
                                 ("const Mesh<Dimension>&", "mesh")],
                         const = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary2(self,
                         args = [("std::vector<Scalar>&", "faceField"),
                                ("const Mesh<Dimension>&", "mesh")],
                         const = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary3(self,
                         args = [("std::vector<Vector>&", "faceField"),
                                 ("const Mesh<Dimension>&", "mesh")],
                         const = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary4(self,
                         args = [("std::vector<Tensor>&", "faceField"),
                                 ("const Mesh<Dimension>&", "mesh")],
                         const = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary5(self,
                         args = [("std::vector<SymTensor>&", "faceField"),
                                 ("const Mesh<Dimension>&", "mesh")],
                         const = True,
                         name = "enforceBoundary"):
        return "void"

    def enforceBoundary6(self,
                         args = [("std::vector<ThirdRankTensor>&", "faceField"),
                                 ("const Mesh<Dimension>&", "mesh")],
                         const = True,
                         name = "enforceBoundary"):
        return "void"

    def swapFaceValues1(self,
                        args = [("Field<Dimension, std::vector<Scalar>>&", "field"),
                                ("const Mesh<Dimension>&", "mesh")],
                        const = True,
                        name = "swapFaceValues"):
        return "void"

    def swapFaceValues2(self,
                        args = [("Field<Dimension, std::vector<Vector>>&", "field"),
                                ("const Mesh<Dimension>&", "mesh")],
                        const = True,
                        name = "swapFaceValues"):
        return "void"

    def meshGhostNodes(self,
                       const = True):
        return "bool"
