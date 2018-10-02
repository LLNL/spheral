#-------------------------------------------------------------------------------
# GenericHydro base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *

@PYB11template("Dimension")
class GenericHydro(Physics):

    typedefs = """
    typedef %(Dimension)s DIM;
    typedef typename DIM::Scalar Scalar;
    typedef typename DIM::Vector Vector;
    typedef typename DIM::Tensor Tensor;
    typedef typename DIM::SymTensor SymTensor;
    typedef typename DIM::ThirdRankTensor ThirdRankTensor;
    typedef typename Physics<DIM>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               W = "const TableKernel<DIM>&",
               WPi = "const TableKernel<DIM>&",
               Q = "ArtificialViscosity<DIM>&",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool"):
        "GenericHydro constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def dt(dataBase = "const DataBase<DIM>&", 
           state = "const State<DIM>&",
           derivs = "const StateDerivatives<DIM>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
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
    artificialViscosity = PYB11property("ArtificialViscosity<DIM>&", "artificialViscosity", doc="The artificial viscosity object")
    kernel = PYB11property("const TableKernel<DIM>&", "kernel", doc="The interpolation kernel")
    PiKernel = PYB11property("const TableKernel<DIM>&", "PiKernel", doc="The interpolation kernel for the artificial viscosity")
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
