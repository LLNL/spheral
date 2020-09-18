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

    PYB11typedefs = """
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename %(Dimension)s::ConvexHull ConvexHull;
  typedef typename %(Dimension)s::FacetedVolume FacetedVolume;
  typedef uint64_t KeyElement;
  typedef std::tuple<KeyElement, KeyElement, KeyElement> Key;
"""

    #...........................................................................
    # Static attributes
    UNSETID = PYB11readonly(static=True, returnpolicy="copy")
    minFacesPerZone = PYB11readonly(static=True, returnpolicy="copy")
    minEdgesPerZone = PYB11readonly(static=True, returnpolicy="copy")
    minNodesPerZone = PYB11readonly(static=True, returnpolicy="copy")
    minEdgesPerFace = PYB11readonly(static=True, returnpolicy="copy")
    minNodesPerFace = PYB11readonly(static=True, returnpolicy="copy")

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
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def valid(self):
        "Perform basic mesh validity checks."
        return "std::string"

    @PYB11virtual
    @PYB11const
    def validDomainInfo(self,
                        xmin = "const Vector&",
                        xmax = "const Vector&",
                        checkUniqueSendProc = "const bool"):
        "Check that the internal parallel info is consistent."
        return "std::string"

    #...........................................................................
    # Static methods
    @PYB11static
    def positiveID(self, id = "const int"):
        "Encapsulate the ones complement for signed (oriented) IDs."
        return "int"

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

    @PYB11pycppname("reconstruct")
    def reconstruct1(self,
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

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    @PYB11pycppname("zone")
    def zone1(self,
              nodeList = "const NodeList<%(Dimension)s>&",
              i = "const unsigned"):
        "We also provide the ability to extract the zone corresponding to the given node in a NodeList."
        return "const Mesh<%(Dimension)s>::Zone&"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    @PYB11pycppname("zone")
    def zone2(self,
              nodeListi = "const unsigned",
              i = "const unsigned"):
        "We also provide the ability to extract the zone corresponding to the given node in a NodeList."
        return "const Mesh<%(Dimension)s>::Zone&"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.nodeBegin(), self.nodeEnd()); }, py::keep_alive<0, 1>()")
    def nodes(self):
        "Return the set of Mesh::Node"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.edgeBegin(), self.edgeEnd()); }, py::keep_alive<0, 1>()")
    def edges(self):
        "Return the set of Mesh::Edge"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.faceBegin(), self.faceEnd()); }, py::keep_alive<0, 1>()")
    def faces(self):
        "Return the set of Mesh::Face"
        return "py::list"

    @PYB11const
    @PYB11implementation("[](const Mesh<%(Dimension)s>& self) { return py::make_iterator(self.zoneBegin(), self.zoneEnd()); }, py::keep_alive<0, 1>()")
    def zones(self):
        "Return the set of Mesh::Zone"
        return "py::list"

    @PYB11const
    def offset(self, nodeList = "const NodeList<%(Dimension)s>&"):
        "Extract the zone offset for the given NodeList."
        return "unsigned"

    @PYB11const
    @PYB11pycppname("offset")
    def offset1(self, nodeListi = "const unsigned"):
        "Extract the zone offset for the given NodeList."
        return "unsigned"

    @PYB11const
    def lookupNodeListID(self,
                         zoneID = "const unsigned",
                         nodeListi = "unsigned&",
                         i = "unsigned&"):
        "Look up the (nodeListID, nodeID) corresponding to the given zoneID."
        return "void"

    def generateDomainInfo(self):
        "Compute the communicated mesh structures."
        return "void"

    def generateParallelRind(self):
        """Generate a parallel rind of cells around each domain representing a one zone
thick set of zones shared with the neighboring processors.
Note we do not recompute the shared elements (nodes & faces) as part of this
procedure, so following this operation those shared elements are no longer
on the surface of the local mesh!"""
        return "void"

    @PYB11pycppname("generateParallelRind")
    def generateParallelRind1(self,
                              generators = "std::vector<Vector>&",
                              Hs = "std::vector<SymTensor>&"):
        "This version also exchanges the generators for the rind cells."
        return "void"

    @PYB11const
    def globalMeshNodeIDs(self):
        "Compute unique global IDs for each node."
        return "std::vector<unsigned>"

    @PYB11const
    def globalMeshFaceIDs(self, globalNodeIDs = "const std::vector<unsigned>&"):
        "Compute unique global IDs for each face."
        return "std::vector<unsigned>"

    def storeNodeListOffsets(self,
                             nodeListPtrs = "const std::vector<NodeList<%(Dimension)s>*>&",
                             offsets = "const std::vector<unsigned>&"):
        "Store the given offsets for a set of NodeLists."
        return "void"

    @PYB11const
    def boundingSurface(self):
        "Compute the bounding surface of the mesh."
        return "FacetedVolume"

    @PYB11const
    def boundingBox(self,
                    xmin = "Vector&",
                    xmax = "Vector&"):
        "Compute the bounding box of the mesh."
        return "void"

    #...........................................................................
    # Properties
    numNodes = PYB11property("unsigned")
    numEdges = PYB11property("unsigned")
    numFaces = PYB11property("unsigned")
    numZones = PYB11property("unsigned")
    neighborDomains = PYB11property("const std::vector<unsigned>&", returnpolicy="reference_internal")
    sharedNodes = PYB11property("const std::vector<std::vector<unsigned>>&", returnpolicy="reference_internal")
    sharedFaces = PYB11property("const std::vector<std::vector<unsigned>>&", returnpolicy="reference_internal")
    minimumScale = PYB11property(doc="Compute the minimum scale (distance between nodes).")

    #---------------------------------------------------------------------------
    # Mesh::Node
    #---------------------------------------------------------------------------
    class Node:
        "Mesh::Node -- a node (or vertex) of a Mesh"

        PYB11typedefs = """
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
        zoneIDs = PYB11property("const std::vector<unsigned>&", returnpolicy="reference_internal")

    #---------------------------------------------------------------------------
    # Mesh::Edge
    #---------------------------------------------------------------------------
    class Edge:
        "Mesh::Edge -- an edge (segment connecting two Nodes) of a Mesh"

        PYB11typedefs = """
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

        PYB11typedefs = """
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

        PYB11typedefs = """
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
