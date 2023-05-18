#-------------------------------------------------------------------------------
# SVPHFactedHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SpheralCommon import *

from GenericHydro import *
from RestartMethods import *

@PYB11template("Dimension")
class SVPHFacetedHydroBase(GenericHydro):
    "SVPHFacetedHydroBase -- The fluid SVPH faceted hydro algorithm"

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename %(Dimension)s::FourthRankTensor FourthRankTensor;
    typedef typename %(Dimension)s::FifthRankTensor FifthRankTensor;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    def pyinit(self,
               smoothingScaleMethod = "const SmoothingScaleBase<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               Q = "ArtificialViscosity<%(Dimension)s>&",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               XSVPH = "const bool",
               linearConsistent = "const bool",
               generateVoid = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               fcentroidal = "const Scalar",
               fcellPressure = "const Scalar",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SVPH Constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual
    @PYB11const
    def dt(self,
           dataBase = "const DataBase<%(Dimension)s>&",
           state = "const State<%(Dimension)s>&",
           derivs = "const StateDerivatives<%(Dimension)s>&",
           currentTime = "const Scalar"):
        "Vote on a time step."
        return "TimeStepType"

    @PYB11virtual 
    def registerState(self,
                      dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(self,
                            dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def initialize(self,
                   time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "void"
                          
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        """Evaluate the derivatives for the principle hydro variables:
mass density, velocity, and specific thermal energy."""
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar",
                 dt = "const Scalar",
                 dataBase = "DataBase<%(Dimension)s>&",
                 state = "State<%(Dimension)s>&",
                 derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the hydro at the completion of an integration step."
        return "void"
                  
    @PYB11virtual
    @PYB11const
    def requireConnectivity(self):
        "This algorithm does not use node->node connectivity."
        return "bool"

    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    #...........................................................................
    # Properties
    kernel = PYB11property(doc="The interpolation kernel")
    densityUpdate = PYB11property("MassDensityType", "densityUpdate", "densityUpdate",
                                  doc="Flag to choose whether we want to sum for density, or integrate the continuity equation.")
    HEvolution = PYB11property("HEvolutionType", "HEvolution", "HEvolution",
                               doc="Flag to select how we want to evolve the H tensor.")
    compatibleEnergyEvolution = PYB11property("bool", "compatibleEnergyEvolution", "compatibleEnergyEvolution",
                                              doc="Flag to determine if we're using the total energy conserving compatible energy evolution scheme.")
    XSVPH = PYB11property("bool", "XSVPH", "XSVPH",
                          doc="Flag to determine if we're using the XSVPH algorithm.")
    linearConsistent = PYB11property("bool", "linearConsistent", "linearConsistent",
                                     doc="Flag to select whether or not to use the linear corrections.")
    generateVoid = PYB11property("bool", "generateVoid", "generateVoid",
                                 doc="Flag to select whether or not to generate void points in the tessellation.")
    fcentroidal = PYB11property("Scalar", "fcentroidal", "fcentroidal",
                                doc="Fraction of centroidal motion to apply each step.")
    fcellPressure = PYB11property("Scalar", "fcellPressure", "fcellPressure",
                                  doc="Fraction of the pressure to take from local cell.")
    xmin = PYB11property("const Vector&", "xmin", "xmin",
                         doc="Optionally we can provide a bounding box for use generating the mesh.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",
                         doc="Optionally we can provide a bounding box for use generating the mesh.")
    smoothingScaleMethod = PYB11property("const SmoothingScaleBase<%(Dimension)s>&", "smoothingScaleMethod",
                                         doc="The object defining how we evolve smoothing scales.")
    mesh = PYB11property("const Mesh<%(Dimension)s>&", "mesh",
                         doc="The tessellation")


    timeStepMask = PYB11property("const FieldList<%(Dimension)s, int>&", "timeStepMask", returnpolicy="reference_internal")
    pressure = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "pressure", returnpolicy="reference_internal")
    cellPressure = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "cellPressure", returnpolicy="reference_internal")
    soundSpeed = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "soundSpeed", returnpolicy="reference_internal")
    volume = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "volume", returnpolicy="reference_internal")
    specificThermalEnergy0 = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "specificThermalEnergy0", returnpolicy="reference_internal")
    Hideal = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "Hideal", returnpolicy="reference_internal")
    maxViscousPressure = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "maxViscousPressure", returnpolicy="reference_internal")
    massDensitySum = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "massDensitySum", returnpolicy="reference_internal")
    weightedNeighborSum = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "weightedNeighborSum", returnpolicy="reference_internal")
    massSecondMoment = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "massSecondMoment", returnpolicy="reference_internal")
    XSVPHDeltaV = PYB11property("const FieldList<%(Dimension)s, Vector>&", "XSVPHDeltaV", returnpolicy="reference_internal")
    DxDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DxDt", returnpolicy="reference_internal")
    DvDt = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DvDt", returnpolicy="reference_internal")
    DmassDensityDt = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DmassDensityDt", returnpolicy="reference_internal")
    DspecificThermalEnergyDt = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DspecificThermalEnergyDt", returnpolicy="reference_internal")
    DHDt = PYB11property("const FieldList<%(Dimension)s, SymTensor>&", "DHDt", returnpolicy="reference_internal")
    DvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "DvDx", returnpolicy="reference_internal")
    internalDvDx = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "internalDvDx", returnpolicy="reference_internal")
    faceForce = PYB11property("const FieldList<%(Dimension)s, std::vector<Vector> >&", "faceForce", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SVPHFacetedHydroBase)

