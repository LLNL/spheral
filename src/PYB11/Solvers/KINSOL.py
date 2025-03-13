#-------------------------------------------------------------------------------
# KINSOL
#-------------------------------------------------------------------------------
from PYB11Generator import *

class KINSOL:
    "Wrapper around the Sundials KINSOL nonlinear solver"

    #...........................................................................
    # Constructors
    def pyinit(self):
        "KINSOL constructor"
        return

    #...........................................................................
    # Methods
    @PYB11virtual
    def solve(self,
              func = "SolverFunction&",
              x = "std::vector<double>&"):
        """Solve the system of equations represented by function 'func' and initial guess 'x'.
Solution vector is returned in 'x'.
Return value indicates success (True) or failure (False)"""
        return "bool"

    #...........................................................................
    # Properties
    globalstrategy = PYB11property("int", "globalstrategy")
    fnormtol = PYB11property("double", "fnormtol")
    scsteptol = PYB11property("double", "scsteptol")
