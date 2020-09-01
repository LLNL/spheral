#-------------------------------------------------------------------------------
# ArtificialViscosity base class
#-------------------------------------------------------------------------------
from PYB11Generator import *
from ArtificialViscosityAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralArtificialViscosity")
class ArtificialViscosity:

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef typename %(Dimension)s::Tensor Tensor;
    typedef typename %(Dimension)s::SymTensor SymTensor;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               Clinear = "const Scalar",
               Cquadratic = "const Scalar",
               QcorrectionOrder = ("const RKOrder", "RKOrder::LinearOrder")):
        "ArtificialViscosity constructor"

    #...........................................................................
    # Methods
    @PYB11const
    def curlVelocityMagnitude(self, DvDx="const Tensor&"):
        "Calculate the curl of the velocity given the stress tensor."
        return "Scalar"

    @PYB11const
    def calculateLimiter(self,
                         vi = "const Vector&",
                         vj = "const Vector&",
                         ci = "const Scalar",
                         cj = "const Scalar",
                         hi = "const Scalar",
                         hj = "const Scalar",
                         nodeListID = "const int",
                         nodeID = "const int"):
        "Method to return the limiter magnitude for the given node."
        return "Tensor"

    @PYB11const
    def shockDirection(self,
                       ci = "const Scalar",
                       hi = "const Scalar",
                       nodeListID = "const int",
                       nodeID = "const int"):
        "Helper for the limiter, calculate the unit grad div v term for the given node"
        return "Vector"

    @PYB11const
    def sigmaWeighting(self, r="const Vector&"):
        "Helper method to calculate the weighting based on the given position for use in the sigma calculation."
        return "Vector"

    @PYB11const
    def sigmaij(self,
                rji = "const Vector&",
                rjiUnit = "const Vector&",
                vji = "const Vector&",
                hi2 = "const Scalar&",
                nodeListID = "const int",
                nodeID = "const int"):
        "Figure out the total stress-strain tensor for a given node pair based on the stored average value and the given (position, velocity) pair."
        return "Tensor"

    #...........................................................................
    # Properties
    Cl = PYB11property("Scalar", "Cl", "Cl",
                       doc="The linear coefficient")
    Cq = PYB11property("Scalar", "Cq", "Cq",
                       doc="The quadratic coefficient")
    QcorrectionOrder = PYB11property("RKOrder", "QcorrectionOrder", "QcorrectionOrder",
                                     doc="The RK correction order used for computing gradients in the viscosity")
    balsaraShearCorrection = PYB11property("bool", "balsaraShearCorrection", "balsaraShearCorrection",
                                           doc="Toggle whether to use the Balsara suppression for shear flows")
    ClMultiplier = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "ClMultiplier",
                                 doc="Correction multiplier for the linear term")
    CqMultiplier = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "CqMultiplier",
                                 doc="Correction multiplier for the quadratic term")
    shearCorrection = PYB11property("const FieldList<%(Dimension)s, Scalar>&", "shearCorrection",
                                    doc="Correction multiplier for Balsara shear suppression")
    sigma = PYB11property("const FieldList<%(Dimension)s, Tensor>&", "sigma",
                          doc="Access the internally computed estimate of sigma: sig^ab = partial v^a / partial x^b")
    gradDivVelocity = PYB11property("const FieldList<%(Dimension)s, Vector>&", "gradDivVelocity",
                                    doc="Access the internally computed estimate of the velocity gradient and grad div velocity")
    limiter = PYB11property("bool", "limiter", "limiter",
                            doc="Toggle whether to apply the del^2 velocity limiter")
    epsilon2 = PYB11property("Scalar", "epsilon2", "epsilon2",
                             doc="Safety factor in denominator for Q")
    negligibleSoundSpeed = PYB11property("Scalar", "negligibleSoundSpeed", "negligibleSoundSpeed",
                                         doc="The negligible sound speed parameter for use in the limiter")
    csMultiplier = PYB11property("Scalar", "csMultiplier", "csMultiplier",
                                 doc="The multiplier for sound speed in the limiter")
    energyMultiplier = PYB11property("Scalar", "energyMultiplier", "energyMultiplier",
                                     doc="The multiplier for energy in the limiter.")

    # This one is a protected property!
    @PYB11const
    @PYB11ignore
    @PYB11protected
    @PYB11cppname("calculateSigma")
    def getcalculateSigma(self):
        return "bool"

    @PYB11ignore
    @PYB11protected
    @PYB11cppname("calculateSigma")
    def setcalculateSigma(self, val="bool"):
        return "void"

    calculateSigma = property(getcalculateSigma, setcalculateSigma, doc="Toggle if sigma should be computed")
    
#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(ArtificialViscosityAbstractMethods, ArtificialViscosity, pure_virtual=True)
PYB11inject(RestartMethods, ArtificialViscosity)
