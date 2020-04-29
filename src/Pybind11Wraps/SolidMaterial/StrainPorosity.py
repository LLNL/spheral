#-------------------------------------------------------------------------------
# StrainPorosity
#-------------------------------------------------------------------------------
from PYB11Generator import *
from Physics import *
from PhysicsAbstractMethods import *
from RestartMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralSolidMaterial")
class StrainPorosity(Physics):
    """An implementation of strain-alpha porosity model described in two papers:

Wunnemann, Collins, & Melosh, Icarus, 180, 514-527 (2006)
"A strain-based porosity model for use in hydrocode simulations of impacts
 and implications for transient crater growth in porous targets"

Collins, G. S., Melosh, H. J., & Wunnemann, K. (2011).
"Improvements to the epsilon-alpha porous compaction model for simulating impacts into
high-porosity solar system objects.
International Journal of Impact Engineering, 38(6), 434-439.
http://doi.org/10.1016/j.ijimpeng.2010.10.013

This model assumes you will provide a solid EOS which will be modified.
The underlying actualy solid EOS should provide the reference density, which
will be treated here as the compacted true solid reference density.

Note this model introduces a new state variable, the distention (alpha), which
the pressure now depends on.  This implies our usual definition of P(rho, eps)
now becomes P(rho, eps, alpha).  Our EOS interface does not recognize this
this parameter, so we store alpha locally and only allow Field updates of the
pressure (forbidding the single value P lookup the EOS usually allows).

StrainPorosity is the physics module which time evolves the distention 
parameter (alpha) and gives it to the PorousEquationOfState."""


    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef typename %(Dimension)s::Vector Vector;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
    typedef typename Physics<%(Dimension)s>::TimeStepType TimeStepType;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               porousEOS = "PorousEquationOfState<%(Dimension)s>&",
               porousStrength = "PorousStrengthModel<%(Dimension)s>&",
               nodeList = "const NodeList<%(Dimension)s>&",
               phi0 = "const double",
               epsE = "const double",
               epsX = "const double",
               kappa = "const double",
               gammaS0 = "const double",
               cS0 = "const double",
               c0 = "const double"):
        """Constructor parameters:
        porousEOS:       Porous EOS we're going to modify
        porousStrength:  Porous strength model we're going to modify
        nodeList:        The NodeList we're going apply to
        phi0:            Initial porosity (single value)
        epsE:            Elastic-plastic transition strain
        epsX:            Threshold strain between compaction regimes
        kappa:           Compaction rate
        gammaS0:         Reference gamma at full density
        cS0:             Reference sound speed at full density
        c0:              Reference sound speed at initial porosity"""

    def pyinit1(self,
                porousEOS = "PorousEquationOfState<%(Dimension)s>&",
                porousStrength = "PorousStrengthModel<%(Dimension)s>&",
                nodeList = "const NodeList<%(Dimension)s>&",
                phi0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&",
                epsE = "const double",
                epsX = "const double",
                kappa = "const double",
                gammaS0 = "const double",
                cS0 = "const double",
                c0 = "const Field<%(Dimension)s, %(Dimension)s::Scalar>&"):
        """Constructor parameters:
        porousEOS:       Porous EOS we're going to modify
        porousStrength:  Porous strength model we're going to modify
        nodeList:        The NodeList we're going apply to
        phi0:            Initial porosity (field of values)
        epsE:            Elastic-plastic transition strain
        epsX:            Threshold strain between compaction regimes
        kappa:           Compaction rate
        gammaS0:         Reference gamma at full density
        cS0:             Reference sound speed at full density
        c0:              Reference sound speed at initial porosity"""

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def initializeProblemStartup(self, 
                                 dataBase = "DataBase<%(Dimension)s>&"):
        "Do any required one-time initializations on problem start up."
        return "void"

    #...........................................................................
    # Properties
    epsE = PYB11property(doc="Elastic-plastic transition strain")
    epsX = PYB11property(doc="Threshold strain between compaction regimes")
    kappa = PYB11property(doc="Compaction rate")
    gammaS0 = PYB11property(doc="Reference gamma at full density")
    cS0 = PYB11property(doc="Reference sound speed at full density")
    porousEOS = PYB11property("const PorousEquationOfState<%(Dimension)s>&", returnpolicy="reference_internal")
    porousStrength = PYB11property("PorousStrengthModel<%(Dimension)s>&", returnpolicy="reference_internal")
    nodeList = PYB11property("const NodeList<%(Dimension)s>&", returnpolicy="reference_internal")
    c0 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    alpha0 = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    alpha = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    DalphaDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    strain = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")
    DstrainDt = PYB11property("const Field<%(Dimension)s, Scalar>&", returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject methods
#-------------------------------------------------------------------------------
PYB11inject(PhysicsAbstractMethods, StrainPorosity, virtual=True, pure_virtual=False)
PYB11inject(RestartMethods, StrainPorosity, virtual=True, pure_virtual=False)
