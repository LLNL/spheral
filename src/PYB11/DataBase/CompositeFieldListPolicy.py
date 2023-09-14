from PYB11Generator import *
from FieldListUpdatePolicyBase import *

@PYB11module("DataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension", "ValueType")
class CompositeFieldListPolicy(FieldListUpdatePolicyBase):
    """CompositeFieldListPolicy -- An implementation of UpdatePolicyBase which
consists of a collection of individual Field policies that should match
the Fields in a FieldList."""

    PYB11typedefs = """
    using KeyType = typename IncrementFieldList<%(Dimension)s, %(ValueType)s>::KeyType;
"""

    #...........................................................................
    # Constructors
    def pyinit0(self):
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

    @PYB11virtual
    def updateAsIncrement(self,
                          key = "const KeyType&",
                          state = "State<%(Dimension)s>&",
                          derivs = "StateDerivatives<%(Dimension)s>&",
                          multiplier = "const double",
                          t = "const double",
                          dt = "const double"):
        """An alternate method to be called when you want to specify that the derivative information
should be assumed to not necessarily be properly time-centered, and therefore you should 
only use time advancement ideas, no "replace" or more sophisticated approaches.
Default to just calling the generic method."""
        return "void"

    def push_back(self,
                  policy = "UpdatePolicyBase<%(Dimension)s>*"):
        "Add a new UpdatePolicy"
        return "void"
