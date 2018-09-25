from PYB11Generator import *

#-------------------------------------------------------------------------------
# NodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
class NodeList:
    "Spheral NodeList base class in %(Dimension)s"

    typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef Field<%(Dimension)s, Vector> VectorField;
    typedef Field<%(Dimension)s, Tensor> TensorField;
    typedef Field<%(Dimension)s, SymTensor> SymTensorField;
"""

    def pyinit(self,
               name = "std::string",
               numInternal = ("unsigned", "0"),
               numGhost = ("unsigned", "0"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("unsigned", "500")):
        "Constructor for NodeList base class."
        return

    @PYB11const
    def mass(self):
        "The mass field"
        return "const ScalarField&"

    @PYB11const
    def positions(self):
        "The position field"
        return "const VectorField&"

    @PYB11const
    def velocity(self):
        "The velocity field"
        return "const VectorField&"

    @PYB11const
    def Hfield(self):
        "The H tensor field"
        return "const SymTensorField&"

    @PYB11const
    def work(self):
        "The CPU work field"
        return "ScalarField&"

    @PYB11pycppname("mass")
    def setmass(self, newValue="const ScalarField&"):
        "Set the mass field"
        return "void"

    @PYB11pycppname("positions")
    def setpositions(self, newValue="const VectorField&"):
        "Set the position field"
        return "void"

    @PYB11pycppname("velocity")
    def setvelocity(self, newValue="const VectorField&"):
        "Set the velocity field"
        return "void"

    @PYB11pycppname("Hfield")
    def setHfield(self, newValue="const SymTensorField&"):
        "Set the H tensor field"
        return "void"

    @PYB11pycppname("work")
    def setwork(self, newValue="const ScalarField&"):
        "Set the CPU work field"
        return "void"

    @PYB11const
    def Hinverse(self, result="SymTensorField&"):
        "Compute the inverse H field"
        return "void"

    @PYB11const
    def haveField(self, field="const FieldBase<%(Dimension)s>&"):
        "Test if the given field is defined on this NodeList"
        return "bool"

    @PYB11const
    def nodeType(self,
                 ID = "int"):
        "Return the classification of the given node"
        return "NodeType"

    @PYB11const
    def neighbor(self):
        "Neighbor object associated with this NodeList"
        return "Neighbor<%(Dimension)s>&"

    def registerNeighbor(self, neighbor="Neighbor<%(Dimension)s>&"):
        "Associate a Neighbor object with this NodeList"
        return "void"

    def unregisterNeighbor(self):
        "Break the relation of this NodeList with it's Neighbor object"
        return "void"

    # Virtual methods
    @PYB11virtual
    def deleteNodes(self, nodeIDs="const std::vector<int>&"):
        "Delete the indicated nodes from this NodeList"
        return "void"

    @PYB11virtual
    def reorderNodes(self, newOrdering="const std::vector<int>&"):
        "Reorder the nodes to the given mapping"
        return "void"

    @PYB11virtual
    @PYB11const
    def label(self):
        "Label for restart files"
        return "std::string"

    @PYB11virtual
    @PYB11const
    def dumpState(self, file="FileIO&", pathName="const std::string&"):
        "Serialize under the given path in a FileIO object"
        return "void"

    @PYB11virtual
    def restoreState(self, file="const FileIO&", pathName="const std::string&"):
        "Restore state from the given path in a FileIO object"
        return "void"

    # Comparison
    def __eq__(self):
        "Equivalence test with another NodeList"

    def __ne__(self):
        "Inequivalence test with another NodeList"

    # Methods used for properties
    @PYB11ignore
    @PYB11cppname("name")
    @PYB11const
    def getname(self):
        return "std::string"

    @PYB11ignore
    @PYB11cppname("numNodes")
    @PYB11const
    def getnumNodes(self):
        return "unsigned"

    @PYB11ignore
    @PYB11cppname("numInternalNodes")
    @PYB11const
    def getnumInternalNodes(self):
        return "unsigned"

    @PYB11ignore
    @PYB11cppname("numGhostNodes")
    @PYB11const
    def getnumGhostNodes(self):
        return "unsigned"

    @PYB11ignore
    @PYB11cppname("numFields")
    @PYB11const
    def getnumFields(self):
        "The number of fields registered on this NodeList"
        return "int"

    @PYB11const
    @PYB11cppname("firstGhostNode")
    def getfirstGhostNode(self):
        "Index of the first ghost node on this NodeList"
        return "unsigned"

    @PYB11ignore
    @PYB11cppname("nodesPerSmoothingScale")
    @PYB11const
    def getnodesPerSmoothingScale(self):
        "Return the target number of nodes per smoothing scale"
        return "Scalar"

    @PYB11cppname("nodesPerSmoothingScale")
    @PYB11ignore
    def setnodesPerSmoothingScale(self, val="const Scalar"):
        "Set the target number of nodes per smoothing scale"
        return "void"

    @PYB11ignore
    @PYB11cppname("maxNumNeighbors")
    @PYB11const
    def getmaxNumNeighbors(self):
        return "unsigned"

    @PYB11ignore
    @PYB11cppname("maxNumNeighbors")
    def setmaxNumNeighbors(self, val="unsigned"):
        return "unsigned"

    @PYB11ignore
    @PYB11cppname("hmin")
    @PYB11const
    def gethmin(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("hmin")
    def sethmin(self, val="double"):
        return "void"

    @PYB11ignore
    @PYB11cppname("hmax")
    @PYB11const
    def gethmax(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("hmax")
    def sethmax(self, val="double"):
        return "void"

    @PYB11ignore
    @PYB11cppname("hminratio")
    @PYB11const
    def gethminratio(self):
        return "double"

    @PYB11ignore
    @PYB11cppname("hminratio")
    def sethminratio(self, val="double"):
        return "void"

    # Properties
    name = property(getname, doc="Name of the NodeList")
    numNodes = property(getnumNodes, doc="Total number of nodes in this NodeList")
    numInternalNodes = property(getnumInternalNodes, doc="Number of internal nodes in this NodeList")
    numGhostNodes = property(getnumGhostNodes, doc="Number of ghost nodes in this NodeListb")
    numFields = property(getnumFields, doc="Number of fields defined on this NodeList")
    firstGhostNode = property(getfirstGhostNode, doc="Index of the first ghost node on this NodeList")
    nodesPerSmoothingScale = property(getnodesPerSmoothingScale, setnodesPerSmoothingScale,
                                      doc="The target number of nodes per smoothing scale")
    maxNumNeighbors = property(getmaxNumNeighbors, setmaxNumNeighbors, doc="The maximum number of neighbors per node allowed for this NodeList")
    hmin = property(gethmin, sethmin, doc="Minimum allowed smoothing scale")
    hmax = property(gethmax, sethmax, doc="Maximum allowed smoothing scale")
    hminratio = property(gethminratio, sethminratio, doc="Minimum allowed ratio of min/max smoothing scale eigenvalues on each node")
