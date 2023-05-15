#-------------------------------------------------------------------------------
# MorrisMonaghanReducingViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class MorrisMonaghanReducingViscosity(Physics):
    """A simple form of the reducing artificial viscosity from Morris & Monaghan.
Computes a correction scaling for the typical Monaghan & Gingold visc.
References:
  Monaghan, J. J, & Gingold, R. A. 1983, J. Comput. Phys., 52, 374
  Monaghan, J. J. 1992, ARA&A, 30, 543
  Morris, J. P., & Monaghan, J. J. 1997, J. Comput. Phys., 136, 41
"""

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
               q = "ArtificialViscosity<%(Dimension)s>&",
               nhQ = ("const Scalar", "5.0"),
               nhL = ("const Scalar", "10.0"),
               aMin = ("const Scalar", "0.1"),
               aMax = ("const Scalar", "2.0")):
        "Morris-Monaghan reducing viscosity evolution constructor"

    #...........................................................................
    # Properties
    DrvAlphaDtQ = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DrvAlphaDtQ")
    DrvAlphaDtL = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "DrvAlphaDtL")
    nhQ = PYB11property("Scalar", "nhQ", "nhQ")
    nhL = PYB11property("Scalar", "nhL", "nhL")
    aMin = PYB11property("Scalar", "aMin", "aMin")
    aMax = PYB11property("Scalar", "aMax", "aMax")

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, MorrisMonaghanReducingViscosity, pure_virtual=False, virtual=True)
PYB11inject(RestartMethods, MorrisMonaghanReducingViscosity)
