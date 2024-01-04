from PYB11Generator import *
from UpdatePolicyBase import *

@PYB11module("SpheralDataBase")
@PYB11holder("std::shared_ptr")
@PYB11template("Dimension")
class FieldUpdatePolicy(UpdatePolicyBase):
    "FieldUpdatePolicy -- Base/interface class for the policies on how Field state variables are to be updated."

    #...........................................................................
    # Virtual methods
    @PYB11virtual
    @PYB11const
    def clonePerField(self):
        "Returns whether this policy should be cloned for each Field in a FieldList or not"
        return "bool"
