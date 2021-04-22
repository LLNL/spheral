#-------------------------------------------------------------------------------
# SolidFSISPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SolidSPHHydroBase import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralSPH")
class SolidFSISPHHydroBase(SolidSPHHydroBase):
    "SolidFSISPHHydroBase -- SolidSPHHydro modified for large density discontinuities"

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
               filter = "const double",
               cfl = "const double",
               surfaceForceCoefficient = "const double",
               densityStabilizationCoefficient = "const double",
               densityDiffusionCoefficient = "const double",
               specificThermalEnergyDiffusionCoefficient = "const double",
               sumDensityNodeLists = "std::vector<int>",
               useVelocityMagnitudeForDt = "const bool",
               compatibleEnergyEvolution = "const bool",
               evolveTotalEnergy = "const bool",
               gradhCorrection = "const bool",
               XSPH = "const bool",
               correctVelocityGradient = "const bool",
               densityUpdate = "const MassDensityType",
               HUpdate = "const HEvolutionType",
               epsTensile = "const double",
               nTensile = "const double",
               damageRelieveRubble = "const bool",
               negativePressureInDamage = "const bool",
               strengthInDamage = "const bool",
               xmin = "const Vector&",
               xmax = "const Vector&"):
        "SolidFSISPHHydroBase constructor"

    #...........................................................................
    # Virtual methods
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


    #...........................................................................
    # Properties
    surfaceNormals = PYB11property("const FieldList<%(Dimension)s, Vector>&", "surfaceNormals", returnpolicy="reference_internal")
    surfaceForceCoefficient = PYB11property("double", "surfaceForceCoefficient", "surfaceForceCoefficient",
                           doc="additional force between different materials ala Monaghan 2013.")
    
    densityStabilizationCoefficient = PYB11property("double", "densityStabilizationCoefficient", "densityStabilizationCoefficient", 
                                          doc="coefficient used to adjust velocity gradient to prevent unstable rho.")
    
    densityDiffusionCoefficient = PYB11property("double", "densityDiffusionCoefficient", "densityDiffusionCoefficient", 
                                          doc="coefficient used to diffuse density amongst like nodes.")
    
    specificThermalEnergyDiffusionCoefficient = PYB11property("double", "specificThermalEnergyDiffusionCoefficient", "specificThermalEnergyDiffusionCoefficient", 
                                          doc="coefficient used to diffuse specificThermalEnergy amongst like nodes.")
    
    sumDensityNodeLists = PYB11property("std::vector<int>", "sumDensityNodeLists", "sumDensityNodeLists", 
                                              doc="control if rigorous density sum is applied to individual node lists.")
    
#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SolidFSISPHHydroBase)
