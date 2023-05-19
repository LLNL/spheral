#-------------------------------------------------------------------------------
# SPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from GenericHydro import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralSPH")
@PYB11dynamic_attr
class SPHHydroBase(GenericHydro):

    PYB11typedefs = """
  typedef typename %(Dimension)s::Scalar Scalar;
  typedef typename %(Dimension)s::Vector Vector;
  typedef typename %(Dimension)s::Tensor Tensor;
  typedef typename %(Dimension)s::SymTensor SymTensor;
  typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""
    
    def pyinit(smoothingScaleMethod = "const SmoothingScaleBase<%(Dimension)s>&",
               dataBase = "DataBase<%(Dimension)s>&",
               Q = "ArtificialViscosity<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               WPi = "const TableKernel<%(Dimension)s>&",
               filter = "const double",
               cfl = "const double",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               gradhCorrection = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               sumMassDensityOverAllNodeLists = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               epsTensile = "const double",
               nTensile = "const double",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SPHHydroBase constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(dataBase = "DataBase<%(Dimension)s>&"):
        "Tasks we do once on problem startup."
        return "void"

    @PYB11virtual 
    def registerState(dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "Register the state Hydro expects to use and evolve."
        return "void"

    @PYB11virtual
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Register the derivatives/change fields for updating state."
        return "void"

    @PYB11virtual
    def preStepInitialize(self,
                          dataBase = "const DataBase<%(Dimension)s>&", 
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Optional hook to be called at the beginning of a time step."
        return "void"

    @PYB11virtual
    def initialize(time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "Initialize the Hydro before we start a derivative evaluation."
        return "void"
                       
    @PYB11virtual
    @PYB11const
    def evaluateDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        """Evaluate the derivatives for the principle hydro 
mass density, velocity, and specific thermal energy."""
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<%(Dimension)s>&",
                            state = "const State<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "Finalize the derivatives."
        return "void"

    @PYB11virtual
    def applyGhostBoundaries(state = "State<%(Dimension)s>&",
                             derivs = "StateDerivatives<%(Dimension)s>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def updateVolume(state = "State<%(Dimension)s>&",
                     boundaries = "const bool"):
        """A method to fill in the volume in the State, optionally enforcing
boundary conditions."""
        return "void"

    #...........................................................................
    # Properties
    kernel = PYB11property("const TableKernel<%(Dimension)s>&", "kernel", doc="The interpolation kernel")
    PiKernel = PYB11property("const TableKernel<%(Dimension)s>&", "PiKernel", doc="The interpolation kernel for the artificial viscosity")
    densityUpdate = PYB11property("MassDensityType", "densityUpdate", "densityUpdate",
                                  doc="Flag to choose whether we want to sum for density, or integrate the continuity equation.")
    HEvolution = PYB11property("HEvolutionType", "HEvolution", "HEvolution",
                               doc="Flag to select how we want to evolve the H tensor")
    compatibleEnergyEvolution = PYB11property("bool", "compatibleEnergyEvolution", "compatibleEnergyEvolution",
                                              doc="Flag to determine if we're using the total energy conserving compatible energy evolution scheme.")
    evolveTotalEnergy = PYB11property("bool", "evolveTotalEnergy", "evolveTotalEnergy",
                                      doc="Flag controlling if we evolve total or specific energy.")
    gradhCorrection = PYB11property("bool", "gradhCorrection", "gradhCorrection",
                                    doc="Flag to determine if we're using the grad h correction.")
    XSPH = PYB11property("bool", "XSPH", "XSPH",
                         doc="Flag to determine if we're using the XSPH algorithm.")
    correctVelocityGradient = PYB11property("bool", "correctVelocityGradient", "correctVelocityGradient",
                                            doc="Flag to determine if we're applying the linear correction for the velocity gradient.")
    sumMassDensityOverAllNodeLists = PYB11property("bool", "sumMassDensityOverAllNodeLists", "sumMassDensityOverAllNodeLists",
                                                   doc="Flag to determine if the sum density definition extends over neighbor NodeLists.")
    filter = PYB11property("double", "filter", "filter", doc="Fraction of position filtering to apply.")
    epsilonTensile = PYB11property("double", "epsilonTensile", "epsilonTensile",
                                   doc="Parameters for the tensile correction force at small scales.")
    nTensile = PYB11property("double", "nTensile", "nTensile",
                             doc="Parameters for the tensile correction force at small scales.")
    xmin = PYB11property("const Vector&", "xmin", "xmin",
                         returnpolicy="reference_internal",
                         doc="Optional minimum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    xmax = PYB11property("const Vector&", "xmax", "xmax",
                         returnpolicy="reference_internal",
                         doc="Optional maximum coordinate for bounding box for use generating the mesh for the Voronoi mass density update.")
    smoothingScaleMethod = PYB11property("const SmoothingScaleBase<%(Dimension)s>&", "smoothingScaleMethod",
                                         returnpolicy="reference_internal",
                                         doc="The object defining how we evolve smoothing scales.")

    timeStepMask =                 PYB11property("const FieldList<%(Dimension)s, int>&",      "timeStepMask",         returnpolicy="reference_internal")
    pressure =                     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "pressure",             returnpolicy="reference_internal")
    soundSpeed =                   PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "soundSpeed",           returnpolicy="reference_internal")
    volume =                       PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "volume",               returnpolicy="reference_internal")
    omegaGradh =                   PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "omegaGradh",           returnpolicy="reference_internal")
    specificThermalEnergy0 =       PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "specificThermalEnergy0",returnpolicy="reference_internal")
    entropy =                      PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "entropy",              returnpolicy="reference_internal")
    Hideal =                       PYB11property("const FieldList<%(Dimension)s, SymTensor>&","Hideal",               returnpolicy="reference_internal")
    maxViscousPressure =           PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "maxViscousPressure",   returnpolicy="reference_internal")
    effectiveViscousPressure =     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "effectiveViscousPressure", returnpolicy="reference_internal")
    massDensityCorrection =        PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "massDensityCorrection",returnpolicy="reference_internal")
    viscousWork =                  PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "viscousWork",          returnpolicy="reference_internal")
    massDensitySum =               PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "massDensitySum",       returnpolicy="reference_internal")
    normalization =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "normalization",        returnpolicy="reference_internal")
    weightedNeighborSum =          PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "weightedNeighborSum",  returnpolicy="reference_internal")
    massSecondMoment =             PYB11property("const FieldList<%(Dimension)s, SymTensor>&","massSecondMoment",     returnpolicy="reference_internal")
    XSPHWeightSum =                PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "XSPHWeightSum",        returnpolicy="reference_internal")
    XSPHDeltaV =                   PYB11property("const FieldList<%(Dimension)s, Vector>&",   "XSPHDeltaV",           returnpolicy="reference_internal")
    M =                            PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "M",                    returnpolicy="reference_internal")
    localM =                       PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "localM",               returnpolicy="reference_internal")
    DxDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DxDt",                 returnpolicy="reference_internal")
    DvDt =                         PYB11property("const FieldList<%(Dimension)s, Vector>&",   "DvDt",                 returnpolicy="reference_internal")
    DmassDensityDt =               PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "DmassDensityDt",       returnpolicy="reference_internal")
    DspecificThermalEnergyDt =     PYB11property("const FieldList<%(Dimension)s, Scalar>&",   "DspecificThermalEnergyDt", returnpolicy="reference_internal")
    DHDt =                         PYB11property("const FieldList<%(Dimension)s, SymTensor>&","DHDt",                 returnpolicy="reference_internal")
    DvDx =                         PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "DvDx",                 returnpolicy="reference_internal")
    internalDvDx =                 PYB11property("const FieldList<%(Dimension)s, Tensor>&",   "internalDvDx",         returnpolicy="reference_internal")
    pairAccelerations =            PYB11property("const std::vector<Vector>&", "pairAccelerations", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SPHHydroBase)
