from PYB11Generator import *
from FieldListUpdatePolicyBase import *

@PYB11module("SpheralHydro")
@PYB11template("Dimension")
@PYB11template_dict({"ValueType" : "typename %(Dimension)s::Scalar"})
class SoundSpeedPolicy(FieldListUpdatePolicyBase):

    PYB11typedefs = """
    using Scalar = typename %(Dimension)s::Scalar;
    using KeyType = typename SoundSpeedPolicy<%(Dimension)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
        return

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    def update(self,
               key = "const KeyType&",
               state = "State<%(Dimension)s>&",
               derivs = "StateDerivatives<%(Dimension)s>&",
               multiplier = "const double",
               t = "const double",
               dt = "const double"):
        "Update a FieldList assoicated with the given key"
        return "void"
