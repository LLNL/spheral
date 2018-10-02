#-------------------------------------------------------------------------------
# CullenDehnenViscosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
class CullenDehnenViscosity(Physics):
    """An implementation of the Cullen Dehnen Viscosity, with Hopkins alterations
Computes a correction for the scaling.
References:
Cullen, L., & Dehnen, W. 2010, MNRAS, 408, 669
Hopkins arXiv:1409.7395
"""

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
               q = "ArtificialViscosity<DIM>&",
               W = "const TableKernel<DIM>&",
               alphMax = ("const Scalar", "2.0"),
               alphMin = ("const Scalar", "0.02"),
               betaC = ("const Scalar", "0.7"),
               betaD = ("const Scalar", "0.05"),
               betaE = ("const Scalar", "1.0"),
               fKern = ("const Scalar", "0.3333333333"),
               boolHopkins = ("const bool", "true")):
        "Cullen-Dehnen viscosity evolution constructor"

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def applyGhostBoundaries(self,
                             state = "State<DIM>&",
                             derivs = "StateDerivatives<DIM>&"):
        "Apply boundary conditions to the physics specific fields."
        return "void"

    @PYB11virtual
    def enforceBoundaries(self,
                          state = "State<DIM>&",
                          derivs = "StateDerivatives<DIM>&"):
        "Enforce boundary conditions for the physics specific fields."
        return "void"

    @PYB11virtual
    def initializeProblemStartup(self,
                                 dataBase = "DataBase<DIM>&"):
        "An optional hook to initialize once when the problem is starting up."
        return "void"

    @PYB11virtual
    def finalize(self,
                 time = "const Scalar", 
                 dt = "const Scalar",
                 dataBase = "DataBase<DIM>&", 
                 state = "State<DIM>&",
                 derivs = "StateDerivatives<DIM>&"):
        "Similarly packages might want a hook to do some post-step finalizations.  Really we should rename this post-step finalize."
        return "void"

    @PYB11virtual
    @PYB11const
    def finalizeDerivatives(self,
                            time = "const Scalar",
                            dt = "const Scalar",
                            dataBase = "const DataBase<DIM>&",
                            state = "const State<DIM>&",
                            derivs = "StateDerivatives<DIM>&"):
        "Provide a hook to be called after all physics packages have had their evaluateDerivatives method called, but before anyone does anything with those derivatives."
        return "void"

    #...........................................................................
    # Properties
    alphMax = PYB11property("Scalar", "alphMax", "alphMax")
    alphMin = PYB11property("Scalar", "alphMin", "alphMin")
    betaE = PYB11property("Scalar", "betaE", "betaE")
    betaD = PYB11property("Scalar", "betaD", "betaD")
    betaC = PYB11property("Scalar", "betaC", "betaC")
    fKern = PYB11property("Scalar", "fKern", "fKern")
    boolHopkins = PYB11property("bool", "boolHopkins", "boolHopkins")
    kernel = PYB11property("const TableKernel<DIM>&", "kernel", returnpolicy="reference_internal")
    PrevDvDt = PYB11property("const FieldList<DIM, Vector>&", "PrevDvDt", returnpolicy="reference_internal")
    PrevDivV = PYB11property("const FieldList<DIM, Scalar>&", "PrevDivV", returnpolicy="reference_internal")
    PrevDivV2 = PYB11property("const FieldList<DIM, Scalar>&", "PrevDivV2", returnpolicy="reference_internal")
    CullAlpha = PYB11property("const FieldList<DIM, Scalar>&", "CullAlpha", returnpolicy="reference_internal")
    CullAlpha2  = PYB11property("const FieldList<DIM, Scalar>&", "CullAlpha2", returnpolicy="reference_internal")
    DalphaDt = PYB11property("const FieldList<DIM, Scalar>&", "DalphaDt", returnpolicy="reference_internal")
    alphaLocal = PYB11property("const FieldList<DIM, Scalar>&", "alphaLocal", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject abstract interface
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, CullenDehnenViscosity, pure_virtual=False, virtual=True)
PYB11inject(RestartMethods, CullenDehnenViscosity)
