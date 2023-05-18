#-------------------------------------------------------------------------------
# HelmholtzEquationOfState
#-------------------------------------------------------------------------------
from PYB11Generator import *
from EquationOfState import *
from EOSAbstractMethods import *

@PYB11template("Dimension")
class HelmholtzEquationOfState(EquationOfState):

    PYB11typedefs = """
    typedef typename %(Dimension)s::Scalar Scalar;
    typedef Field<%(Dimension)s, Scalar> ScalarField;
"""

    #...........................................................................
    # Constructors
    def pyinit(self,
               constants = "const PhysicalConstants&",
               minimumPressure = ("const double", "std::numeric_limits<double>::lowest()"),
               maximumPressure = ("const double",  "std::numeric_limits<double>::max()"),
               minimumTemperature = ("const double", "std::numeric_limits<double>::lowest()"),
               minPressureType = ("const MaterialPressureMinType", "MaterialPressureMinType::PressureFloor"),
               abar0 = ("double", "13.6"),
               zbar0 = ("double", "6.8"),
               externalPressure = ("double", "0.0")):
        "Helmholtz constructor"

    #...........................................................................
    # Properties
    abar = PYB11property(returnpolicy="reference_internal")
    zbar = PYB11property(returnpolicy="reference_internal")
    needUpdate = PYB11property("bool", "getUpdateStatus", "setUpdateStatus")
    
#-------------------------------------------------------------------------------
# Add the virtual interface
#-------------------------------------------------------------------------------
PYB11inject(EOSAbstractMethods, HelmholtzEquationOfState, virtual=True)
