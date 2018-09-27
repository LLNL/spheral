#-------------------------------------------------------------------------------
# StateDerivatives
#-------------------------------------------------------------------------------
from PYB11Generator import *
from DataBase import *

@PYB11template("Dimension")
class DataBase:

    typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename DataBase<%(Dimension)s>::ConnectivityMapPtr ConnectivityMapPtr;
"""

    #...........................................................................
    # Constructors
    def pyinit(self):
        "Default constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def reinitializeNeighbors(self):
        "Optimize all Neighbor objects for the current state."
        return "void"

    @PYB11const
    def updateConnectivityMap(self, computeGhostConnectivity="const bool"):
        "Update the internal connectivity map."
        return "void"

    @PYB11const
    def patchConnectivityMap(self,
                             flags = "const FieldList<%(Dimension)s, int>&",
                             old2new = "const FieldList<%(Dimension)s, int>&"):
        "Update the internal connectivity map."
        return "void"

    @PYB11returnpolicy("reference_internal")
    @PYB11const
    def connectivityMap(self, computeGhostConnectivity=("const bool", "false")):
        "Get the connectivity map, optionally including ghost connectivity"
        return "const ConnectivityMap<%(Dimension)s>&"

    @PYB11const
    def connectivityMapPtr(self, computeGhostConnectivity=("const bool", "false")):
        "Get the connectivity map as a std::shared_ptr, optionally including ghost connectivity"
        return "ConnectivityMapPtr"

    #...........................................................................
    # Properties
    numNodeLists = PYB11property("int", "numNodeLists", doc="Number of NodeLists in DataBase")
    numFluidNodeLists = PYB11property("int", "numFluidNodeLists", doc="Number of FluidNodeLists in DataBase")
    numSolidNodeLists = PYB11property("int", "numSolidNodeLists", doc="Number of SolidNodeLists in DataBase")

    numNodes = PYB11property("int", "numNodes", doc="Number of nodes in all NodeLists in DataBase")
    numInternalNodes = PYB11property("int", "numInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase")
    numGhostNodes = PYB11property("int", "numGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase")

    globalNumNodes = PYB11property("int", "globalNumNodes", doc="Number of nodes in all NodeLists in DataBase across all processors")
    globalNumInternalNodes = PYB11property("int", "globalNumInternalNodes", doc="Number of internal nodes in all NodeLists in DataBase across all processors")
    globalNumGhostNodes = PYB11property("int", "globalNumGhostNodes", doc="Number of ghost nodes in all NodeLists in DataBase across all processors")
