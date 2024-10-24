#-------------------------------------------------------------------------------
# GenericHydro base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
@PYB11module("SpheralPhysics")
class GenericHydro(Physics):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Q = "ArtificialViscosity<%(Dimension)s>&",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool"):
        "GenericHydro constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<%(Dimension)s>&", 
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step (explicit)"
        return "TimeStepType"

    @PYB11virtual
    @PYB11const
    def dtImplicit(dataBase = "const DataBase<%(Dimension)s>&", 
                   state = "const State<%(Dimension)s>&",
                   derivs = "const StateDerivatives<%(Dimension)s>&",
                   currentTime = "const Scalar"):
        "Vote on a time step (implicit)"
        return "TimeStepType"

    #...........................................................................
    # Protected methods
    @PYB11protected
    @PYB11const
    def updateMasterNeighborStats(self, numMaster="int"):
        return "void"

    @PYB11protected
    @PYB11const
    def updateCoarseNeighborStats(self, numCoarse="int"):
        return "void"

    @PYB11protected
    @PYB11const
    def updateRefineNeighborStats(self, numRefine="int"):
        return "void"

    @PYB11protected
    @PYB11const
    def updateActualNeighborStats(self, numActual="int"):
        return "void"

    #...........................................................................
    # Attributes
    artificialViscosity = PYB11property("ArtificialViscosity<%(Dimension)s>&", "artificialViscosity", doc="The artificial viscosity object")
    cfl = PYB11property("Scalar", "cfl", "cfl", doc="The Courant-Friedrichs-Lewy timestep limit multiplier")
    useVelocityMagnitudeForDt = PYB11property("bool", "useVelocityMagnitudeForDt", "useVelocityMagnitudeForDt", doc="Should the pointwise velocity magnitude be used to limit the timestep?")
    minMasterNeighbor = PYB11property("int", "minMasterNeighbor", doc="minimum number of master neighbors found")
    maxMasterNeighbor = PYB11property("int", "maxMasterNeighbor", doc="maximum number of master neighbors found")
    averageMasterNeighbor = PYB11property("double", "averageMasterNeighbor", doc="average number of master neighbors found")
    minCoarseNeighbor = PYB11property("int", "minCoarseNeighbor", doc="minimum number of coarse neighbors found")
    maxCoarseNeighbor = PYB11property("int", "maxCoarseNeighbor", doc="maximum number of coarse neighbors found")
    averageCoarseNeighbor = PYB11property("double", "averageCoarseNeighbor", doc="average number of coarse neighbors found")
    minRefineNeighbor = PYB11property("int", "minRefineNeighbor", doc="minimum number of refine neighbors found")
    maxRefineNeighbor = PYB11property("int", "maxRefineNeighbor", doc="maximum number of refine neighbors found")
    averageRefineNeighbor = PYB11property("double", "averageRefineNeighbor", doc="average number of refine neighbors found")
    minActualNeighbor = PYB11property("int", "minActualNeighbor", doc="minimum number of actual neighbors found")
    maxActualNeighbor = PYB11property("int", "maxActualNeighbor", doc="maximum number of actual neighbors found")
    averageActualNeighbor = PYB11property("double", "averageActualNeighbor", doc="average number of actual neighbors found")
    DTrank = PYB11property("size_t", "DTrank", doc="rank of processor controlling last time step")
    DTNodeList = PYB11property("size_t", "DTNodeList", doc="NodeList index of NodeList controlling last time step")
    DTnode = PYB11property("size_t", "DTnode", doc="Node ID of node controlling last time step")
    DTreason = PYB11property("std::string", "DTreason", doc="short string describing constraint controlling last time step")
