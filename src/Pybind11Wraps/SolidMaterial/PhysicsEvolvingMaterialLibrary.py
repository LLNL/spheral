#-------------------------------------------------------------------------------
# PhysicsEvolvingMaterialLibrary
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from SolidEquationOfState import *
from StrengthModel import *

@PYB11template("Dimension")
class PhysicsEvolvingMaterialLibrary(Physics,
                                     SolidEquationOfState,
                                     StrengthModel):
    """PhysicsEvolvingMaterialLibrary -- A base class to combine the Physics,
 EquationOfState, and StrengthModel for more sophisticated material
 libraries.

This class implements three distinct Spheral interfaces:
  EquationOfState : provides the ordinary Spheral EOS interface
  StrengthModel   : also answers the Strength questions about shear modulus 
                    and yield strength
  Physics         : Geodyn also needs to advance it's own internal state, so
                    this class implements the Physics interface to support
                    that."""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
    typedef typename %(Dimension)s::ThirdRankTensor ThirdRankTensor;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
"""


    def pyinit(self,
               referenceDensity = "const double",
               etamin = "const double",
               etamax = "const double",
               constants = "const PhysicalConstants&",
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double", "std::numeric_limits<double>::max()"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor")):
        "Constructor"


