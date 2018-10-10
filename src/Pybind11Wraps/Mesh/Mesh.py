#-------------------------------------------------------------------------------
# Mesh class
#-------------------------------------------------------------------------------
from PYB11Generator import *

@PYB11template("Dimension")
class Mesh:
    """Mesh -- the generic mesh class we construct around meshless nodes.
1D:  Mesh is a collection of line segments.
2D:  Arbitrary polygons.
3D:  Arbitrary polyhedra."""

    typedefs = """
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename %(Dimension)s::ConvexHull ConvexHull;
  typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
  typedef uint64_t KeyElement;
  typedef boost::tuple<KeyElement, KeyElement, KeyElement> Key;
"""

    #...........................................................................
    # Static attributes
    UNSETID = PYB11readonly(static=True)
    minFacesPerZone = PYB11readonly(static=True)
    minEdgesPerZone = PYB11readonly(static=True)
    minNodesPerZone = PYB11readonly(static=True)
    minEdgesPerFace = PYB11readonly(static=True)
    minNodesPerFace = PYB11readonly(static=True)

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    def pyinit1(self,
                generators = "const std::vector<Vector>&",
                xmin = "const Vector&",
                xmax = "const Vector&"):
        "Construct as a tesselation of points in a box"

    def pyinit2(self,
                generators = "const std::vector<Vector>&",
                boundary = "const FacetedVolume&"):
        "Construct as a tesselation of points in a bounding volume"

    def pyinit3(self,
                nodePositions = "const std::vector<Vector>&",
                edgeNodes = "const std::vector<std::vector<unsigned> >&",
                faceEdges = "const std::vector<std::vector<unsigned> >&",
                zoneFaces = "const std::vector<std::vector<int> >&"):
        "Construct with explicitly specified geometry & topology"

    #...........................................................................
    # Methods
    def clear(self):
        "Clear out any exising data in the mesh."
        return "void"

    def reconstruct(self,
                    generators = "const std::vector<Vector>&",
                    xmin = "const Vector&",
                    xmax = "const Vector&"):
        "Reconstruct from a set of generators in a box."
        return "void"

    def reconstruct(self,
                    generators = "const std::vector<Vector>&",
                    boundary = "const FacetedVolume&"):
        "Reconstruct from a set of generators in a bounding volume."
        return "void"

    def removeZonesByMask(self,
                          mask = "const std::vector<unsigned>&"):
        """Remove zones from the mesh according to a mask:
   mask[i] = 0  ---> remove zone i
   mask[i] = 1  ---> keep zone i"""
        return "void"

    def cleanEdges(self, edgeTol = "const double"):
        "Remove edges below a threshold fraction size."
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def node(self, i="const unsigned"):
        "Return a Mesh::Node by index"
        return "const Mesh<%(Dimension)s>::Node&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def edge(self, i="const unsigned"):
        "Return a Mesh::Edge by index"
        return "const Mesh<%(Dimension)s>::Edge&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def face(self, i="const unsigned"):
        "Return a Mesh::Face by index"
        return "const Mesh<%(Dimension)s>::Face&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def zone(self, i="const unsigned"):
        "Return a Mesh::Zone by index"
        return "const Mesh<%(Dimension)s>::Zone&"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.nodeBegin(), self.nodeEnd()); }, py::keep_alive<0, 1>()")
    def nodes(self):
        "Return a Mesh::Node by index"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.edgeBegin(), self.edgeEnd()); }, py::keep_alive<0, 1>()")
    def edges(self):
        "Return a Mesh::Edge by index"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.faceBegin(), self.faceEnd()); }, py::keep_alive<0, 1>()")
    def faces(self):
        "Return a Mesh::Face by index"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.zoneBegin(), self.zoneEnd()); }, py::keep_alive<0, 1>()")
    def zones(self):
        "Return a Mesh::Zone by index"
        return "py::list"

    #...........................................................................
    # Properties
    numNodes = PYB11property("unsigned")
    numEdges = PYB11property("unsigned")
    numFaces = PYB11property("unsigned")
    numZones = PYB11property("unsigned")

    #---------------------------------------------------------------------------
    # Mesh::Node
    #---------------------------------------------------------------------------
    class Node:
        "Mesh::Node -- a node (or vertex) of a Mesh"

        typedefs = """
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

        def pyinit(self,
                   mesh = "const Mesh<%(Dimension)s>&",
                   ID = "const unsigned",
                   zoneIDs = "const std::vector<unsigned>&"):
            "Constructor"

        ID = PYB11property("unsigned")
        position = PYB11property("Vector")
        zoneIDs = PYB11property("const std::vector<unsigned>&")

    #---------------------------------------------------------------------------
    # Mesh::Edge
    #---------------------------------------------------------------------------
    class Edge:
        "Mesh::Edge -- an edge (segment connecting two Nodes) of a Mesh"

        typedefs = """
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

        def pyinit(self,
                   mesh = "const Mesh<%(Dimension)s>&",
                   ID = "const unsigned",
                   node1ID = "const unsigned",
                   node2ID = "const unsigned"):
            "Constructor"

        ID = PYB11property("unsigned")
        node1ID = PYB11property("unsigned")
        node2ID = PYB11property("unsigned")
        node1 = PYB11property("const Mesh<%(Dimension)s>::Node&", returnpolicy="reference_internal")
        node2 = PYB11property("const Mesh<%(Dimension)s>::Node&", returnpolicy="reference_internal")
        position = PYB11property("Vector")
        length = PYB11property("double")

    #---------------------------------------------------------------------------
    # Mesh::Face
    #---------------------------------------------------------------------------
    class Face:
        "Mesh::Face -- an planar facet of a Zone bounded by a closed loop of Edges"

        typedefs = """
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

        def pyinit(self,
                   mesh = "const Mesh<%(Dimension)s>&",
                   ID = "const unsigned",
                   zone1ID = "const int",
                   zone2ID = "const int",
                   edgeIDs = "const std::vector<unsigned>&"):
            "Constructor"

        #......................................................................
        # Methods
        @PYB11const
        def oppositeZoneID(self, zoneID = "const int"):
            "Return the other zone sharing this face."
            return "int"

        def compare(self,
                    point = "const Vector&",
                    tol = ("const double", "1.0e-8")):
            "Is the given point above, below, or coplanar with the facet?"
            return "int"

        #......................................................................
        # Properties
        ID = PYB11property("unsigned")
        numNodes = PYB11property("unsigned")
        numEdges = PYB11property("unsigned")
        nodeIDs = PYB11property("const std::vector<unsigned>&", returnpolicy="reference_internal")
        edgeIDs = PYB11property("const std::vector<unsigned>&", returnpolicy="reference_internal")
        zone1ID = PYB11property("int")
        zone2ID = PYB11property("int")
        position = PYB11property("Vector")
        area = PYB11property("double")
        unitNormal = PYB11property("Vector")

    #---------------------------------------------------------------------------
    # Mesh::Zone
    #---------------------------------------------------------------------------
    class Zone:
        "Mesh::Zone -- a volume bounded by planar Faces"

        typedefs = """
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::SymTensor SymTensor;
"""

        def pyinit(self,
                   mesh = "const Mesh<%(Dimension)s>&",
                   ID = "const unsigned",
                   faceIDs = "const std::vector<int>&"):
            "Constructor"

        ID = PYB11property("unsigned")
        numNodes = PYB11property("unsigned")
        numEdges = PYB11property("unsigned")
        numFaces = PYB11property("unsigned")
        nodeIDs = PYB11property("const std::vector<unsigned>&", returnpolicy="reference_internal")
        edgeIDs = PYB11property("const std::vector<unsigned>&", returnpolicy="reference_internal")
        faceIDs = PYB11property("const std::vector<int>&", returnpolicy="reference_internal")
        position = PYB11property("Vector")
        volume = PYB11property("double")
        convexHull = PYB11property("typename Mesh<%(Dimension)s>::ConvexHull")
