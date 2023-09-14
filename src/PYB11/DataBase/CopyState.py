from PYB11Generator import *
from FieldUpdatePolicyBase import *

@PYB11module("SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension", "ValueType")
class CopyState(FieldUpdatePolicyBase):
    """CopyState -- An implementation of UpdatePolicyBase appropriate for
copying one state field to another.
Assumes that the copied state is dependent upon the master state, *and* that
that is the only dependency."""

    PYB11typedefs = """
    using KeyType = typename CopyState<%(Dimension)s, %(ValueType)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self,
                masterState = "const std::string&",
                copyState = "const std::string&"):
        return

    #...........................................................................
    # Methods
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

    @PYB11static
    def prefix(self):
        "Prefix for key of derivatives"
        return "const std::string"
