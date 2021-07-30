from PYB11Generator import *
from NodeList import NodeList
from RestartMethods import *

#-------------------------------------------------------------------------------
# NeighborNodeList template
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralNodeList")
@PYB11dynamic_attr
class NeighborNodeList(NodeList):
    "Spheral NeighborNodeList base class in %(Dimension)s, i.e.,  the NodeList for neighbors."

    PYB11typedefs = """
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
               numInternal = ("int", "0"),
               numGhost = ("int", "0"),
               hmin = ("double", "1e-20"),
               hmax = ("double", "1e20"),
               hminratio = ("double", "0.1"),
               nPerh = ("double", "2.01"),
               maxNumNeighbors = ("int", "500")):
        "Constructor for a NeighborNodeList class."
        return

    @PYB11const
    @PYB11returnpolicy("reference_internal")
    def neighbor(self):
        "Neighbor object associated with this NodeList"
        return "Neighbor<%(Dimension)s>&"

    def registerNeighbor(self, neighbor="Neighbor<%(Dimension)s>&"):
        "Associate a Neighbor object with this NodeList"
        return "void"

    def unregisterNeighbor(self):
        "Break the relation of this NodeList with it's Neighbor object"
        return "void"

    #...........................................................................
    # Comparison
    def __eq__(self):
        "Equivalence test with another NeighborNodeList"

    def __ne__(self):
        "Inequivalence test with another NeighborNodeList"

    #...........................................................................
    # Properties
    maxNumNeighbors = PYB11property("unsigned", "maxNumNeighbors", "maxNumNeighbors",
                                    doc="The maximum number of neighbors per node allowed for this NodeList")
    
#-------------------------------------------------------------------------------
# Inject the restart methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, NeighborNodeList)
