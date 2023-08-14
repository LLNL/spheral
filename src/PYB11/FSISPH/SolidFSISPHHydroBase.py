#-------------------------------------------------------------------------------
# SolidFSISPHHydroBase
#-------------------------------------------------------------------------------
from PYB11Generator import *
from SPHHydroBase import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralFSISPH")
class SolidFSISPHHydroBase(SPHHydroBase):
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
               slides = "SlideSurface<%(Dimension)s>&",
               W = "const TableKernel<%(Dimension)s>&",
               filter = "const double",
               cfl = "const double",
               surfaceForceCoefficient = "const double",
               densityStabilizationCoefficient = "const double",
               specificThermalEnergyDiffusionCoefficient = "const double",
               xsphCoefficient = "const double",
               interfaceMethod = "const InterfaceMethod",
               kernelAveragingMethod = "const KernelAveragingMethod",
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

    @PYB11virtual
    def initializeProblemStartup(dataBase = "DataBase<%(Dimension)s>&"):
        "register the surface normals w/ the database"
        return "void"


    @PYB11virtual
    def initialize(time = "const Scalar",
                   dt = "const Scalar",
                   dataBase = "const DataBase<%(Dimension)s>&",
                   state = "State<%(Dimension)s>&",
                   derivs = "StateDerivatives<%(Dimension)s>&"):
        "calculates surface normals"
        return "void"

    
    @PYB11virtual
    def registerState(dataBase = "DataBase<%(Dimension)s>&",
                      state = "State<%(Dimension)s>&"):
        "register the surface normals"
        return "void"


    @PYB11virtual
    def registerDerivatives(dataBase = "DataBase<%(Dimension)s>&",
                            derivs = "StateDerivatives<%(Dimension)s>&"):
        "non-op place filler"
        return "void"



    #...........................................................................
    # Properties
    DepsDx = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DepsDx", returnpolicy="reference_internal")
    DPDx = PYB11property("const FieldList<%(Dimension)s, Vector>&", "DPDx", returnpolicy="reference_internal")
    rawPressure = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "rawPressure", returnpolicy="reference_internal")
    interfaceNormals = PYB11property("const FieldList<%(Dimension)s, Vector>&", "interfaceNormals", returnpolicy="reference_internal")
    interfaceFraction = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "interfaceFraction", returnpolicy="reference_internal")
    interfaceSmoothness = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "interfaceSmoothness", returnpolicy="reference_internal")
    newInterfaceNormals = PYB11property("const FieldList<%(Dimension)s, Vector>&", "newInterfaceNormals", returnpolicy="reference_internal")
    newInterfaceFraction = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "newInterfaceFraction", returnpolicy="reference_internal")
    newInterfaceSmoothness = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "newInterfaceSmoothness", returnpolicy="reference_internal")
    
    slideSurfaces = PYB11property("SlideSurface<%(Dimension)s>&", "slideSurface", doc="The slide surface object")
    
    surfaceForceCoefficient = PYB11property("double", "surfaceForceCoefficient", "surfaceForceCoefficient",doc="additional force between different materials ala Monaghan 2013.")
    densityStabilizationCoefficient = PYB11property("double", "densityStabilizationCoefficient", "densityStabilizationCoefficient",doc="coefficient used to adjust velocity gradient to prevent unstable rho.")
    specificThermalEnergyDiffusionCoefficient = PYB11property("double", "specificThermalEnergyDiffusionCoefficient", "specificThermalEnergyDiffusionCoefficient",doc="coefficient used to diffuse specificThermalEnergy amongst like nodes.")
    xsphCoefficient = PYB11property("double", "xsphCoefficient", "xsphCoefficient",doc="coefficient to dial magnitude of xsph.")
    
    sumDensityNodeLists = PYB11property("std::vector<int>", "sumDensityNodeLists", "sumDensityNodeLists",doc="control if rigorous density sum is applied to individual node lists.")
    
    interfaceMethod = PYB11property("InterfaceMethod", "interfaceMethod", "interfaceMethod",doc="Flag to select how we want construct material interfaces")
    kernelAveragingMethod = PYB11property("KernelAveragingMethod", "kernelAveragingMethod", "kernelAveragingMethod",doc="Flag to select our kernel type")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(RestartMethods, SolidFSISPHHydroBase)
