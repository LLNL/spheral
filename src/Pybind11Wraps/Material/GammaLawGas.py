#-------------------------------------------------------------------------------
# GammaLawGas
#-------------------------------------------------------------------------------
from PYB11Generator import *
from EquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
class GammaLawGas(EquationOfState):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               gamma = "const double",
               mu = "const double",
               constants = "const PhysicalConstants&",
               minimumPressure = ("const double", "-std::numeric_limits<double>::max()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor")):
        "Gamma law gas constructor: gamma=ratio of specific heats, mu=mean molecular weight"

    #...........................................................................
    # Methods
    @PYB11const
    def pressure(self,
                 massDensity = "const Scalar",
                 specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def temperature(self,
                    massDensity = "const Scalar",
                    specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def specificThermalEnergy(self,
                              massDensity = "const Scalar",
                              temperature = "const Scalar"):
        return "Scalar"

    @PYB11const
    def specificHeat(self,
                     massDensity = "const Scalar",
                     temperature = "const Scalar"):
        return "Scalar"

    @PYB11const
    def soundSpeed(self,
                   massDensity = "const Scalar",
                   specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    @PYB11pycppname("gamma")
    def gamma1(self,
               massDensity = "const Scalar",
               specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def bulkModulus(self,
                    massDensity = "const Scalar",
                    specificThermalEnergy = "const Scalar"):
        return "Scalar"

    @PYB11const
    def entropy(self,
                massDensity = "const Scalar",
                specificThermalEnergy = "const Scalar"):
        return "Scalar"

    #...........................................................................
    # Properties
    gamma = PYB11property("double", "gamma", "gamma", doc="gamma: ratio of specific heats")
    mu = PYB11property("double", "molecularWeight", "molecularWeight", doc="mean molecular weight")
    
#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, GammaLawGas, virtual=True)
