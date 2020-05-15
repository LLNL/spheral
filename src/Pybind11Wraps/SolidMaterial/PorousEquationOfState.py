#-------------------------------------------------------------------------------
# PorousEquationOfState
#-------------------------------------------------------------------------------
from PYB11Generator import *
from EquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
@PYB11module("SpheralSolidMaterial")
class PorousEquationOfState(EquationOfState):
    """An implementation of strain-alpha porosity model described in
Wunnemann, Collins, & Melosh, Icarus, 180, 514-527 (2006)
"A strain-based porosity model for use in hydrocode simulations of impacts
 and implications for transient crater growth in porous targets"

This model assumes you will provide a solid EOS which will be modified.
The underlying actualy solid EOS should provide the reference density, which
will be treated here as the compacted true solid reference density.

Note this model introduces a new state variable, the distention (alpha), which
the pressure now depends on.  This implies our usual definition of P(rho, eps)
now becomes P(rho, eps, alpha).  Our EOS interface does not recognize this
this parameter, so we store alpha locally and only allow Field updates of the
pressure (forbidding the single value P lookup the EOS usually allows)."""

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               solidEOS = "const EquationOfState<%(Dimension)s>&"):
        "Porous EOS constructor"

    #...........................................................................
    # Properties
    solidEOS = PYB11property("const EquationOfState<%(Dimension)s>&", returnpolicy="reference_internal")
    alpha = PYB11property("const Field<%(Dimension)s, Scalar>&", "alpha", "alpha",
                          returnpolicy="reference_internal")
    alpha0 = PYB11property("const Field<%(Dimension)s, Scalar>&", "alpha0", "alpha0",
                           returnpolicy="reference_internal")
    c0 = PYB11property("const Field<%(Dimension)s, Scalar>&", "c0", "c0",
                       returnpolicy="reference_internal")

#-------------------------------------------------------------------------------
# Inject the EOS interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, PorousEquationOfState, virtual=True, pure_virtual=False)
