#-------------------------------------------------------------------------------
# PalphaPorosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralPorosity")
class PalphaPorosity(Physics):
    """An implementation of the P-alpha porosity model described in

Jutzi, M., Benz, W., & Michel, P. (2008). Numerical simulations of impacts
involving porous bodies.I. Implementing sub-resolution porosity in a 3D SPH
hydrocode. Icarus, 198(1), 242–255.

This model assumes you will provide a solid EOS which will be modified.
The underlying actualy solid EOS should provide the reference density, which
will be treated here as the compacted true solid reference density.

Note this model introduces a new state variable, the distention (alpha), which
the pressure now depends on.  This implies our usual definition of P(rho, eps)
now becomes P(rho, eps, alpha).  Our EOS interface does not recognize this
this parameter, so we store alpha locally and only allow Field updates of the
pressure (forbidding the single value P lookup the EOS usually allows).

PalphaPorosity is the physics module which time evolves the distention 
parameter (alpha) and gives it to the PorousEquationOfState.
"""


    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               nodeList = "const SolidNodeList<%(Dimension)s>&",
               phi0 = "const double",
               Pe = "const double",
               Pt = "const double",
               Ps = "const double",
               alphae = "const double",
               alphat = "const double",
               n1 = "const double",
               n2 = "const double",
               cS0 = "const double",
               c0 = "const double",
               rho0 = "const double"):
        """Constructor parameters:
        nodeList:        The SolidNodeList we're going apply to
        phi0:            Initial porosity (single value)
        Pe:              Elastic pressure threshold
        Pt:              Transition pressure (Pe <= Pt)
        Ps:              Solid transition pressure (from fit, Pt <= Ps)
        alphae:          alpha at P=Pe
        alphat:          Transition distension
        n1:              Fitted exponent for plastic distention evolution
        n2:              Fitted exponent for plastic distention evolution
        cS0:             Reference sound speed at full density
        c0:              Reference sound speed at initial porosity
        rho0:            Reference solid density"""

    def pyinit1(self,
                nodeList = "const SolidNodeList<%(Dimension)s>&",
                phi0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                Pe = "const double",
                Pt = "const double",
                Ps = "const double",
                alphae = "const double",
                alphat = "const double",
                n1 = "const double",
                n2 = "const double",
                cS0 = "const double",
                c0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                rho0 = "const double"):
        """Constructor parameters:
        nodeList:        The SolidNodeList we're going apply to
        phi0:            Initial porosity (field of values)
        Pe:              Elastic pressure threshold
        Pt:              Transition pressure (Pe <= Pt)
        Ps:              Solid transition pressure (from fit, Pt <= Ps)
        alphae:          alpha at P=Pe
        alphat:          Transition distension
        n1:              Fitted exponent for plastic distention evolution
        n2:              Fitted exponent for plastic distention evolution
        cS0:             Reference sound speed at full density
        c0:              Reference sound speed at initial porosity
        rho0:            Reference solid density"""

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, 
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Do any required one-time initializations on problem start up."
        return "void"

    #...........................................................................
    # Methods
    @PYB11const
    def phi(self):
        "Compute the current porosity from the distention"
        return "ScalarField"

    #...........................................................................
    # Properties
    Pe = PYB11property(doc="Elastic pressure threshold")
    Pt = PYB11property(doc="Transition pressure (Pe <= Pt)")
    Ps = PYB11property(doc="Transition pressure (Pe <= Pt)")
    alphae = PYB11property(doc="Elastic distension threshold")
    alphat = PYB11property(doc="Elastic distension threshold")
    n1 = PYB11property(doc="Fitted exponent for plastic distention evolution")
    n2 = PYB11property(doc="Fitted exponent for plastic distention evolution")
    cS0 = PYB11property(doc="Reference sound speed at full density")
    K0 = PYB11property()
    rho0 = PYB11property(doc="Reference solid density")
    fdt = PYB11property("double", getter="fdt", setter="fdt", doc="The timestep fractional multiplier (0 => no timestep control on alpha)")
    maxAbsDalphaDt = PYB11property(doc="maximum of the last abs(DalphaDt) calculated")
    nodeList = PYB11property("const SolidNodeList<%(Dimension)s>&", returnpolicy="reference_internal")
    c0 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    alpha0 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    alpha = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    DalphaDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    partialPpartialEps = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    partialPpartialRho = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, PalphaPorosity, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, PalphaPorosity, virtual=True, pure_virtual=False)
