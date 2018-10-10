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
