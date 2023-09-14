from PYB11Generator import *
from FieldListUpdatePolicyBase import *

@PYB11template("Dimension", "ValueType")
class CopyFieldList(FieldListUpdatePolicyBase):
    """CopyFieldList -- An implementation of UpdatePolicyBase appropriate for
copying one state field to another.
Assumes that the copied state is dependent upon the master state, *and* that
that is the only dependency."""

    PYB11typedefs = """
    using KeyType = typename IncrementFieldList<%(Dimension)s, %(ValueType)s>::KeyType;
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
