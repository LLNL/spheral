from PYB11Generator import *

#-------------------------------------------------------------------------------
# FieldSpanBase
#-------------------------------------------------------------------------------
@PYB11template("Dimension")
@PYB11module("SpheralField")
class FieldSpanBase:

    @PYB11pure_virtual
    @PYB11const
    def size(self):
        "Number of elements"
        return "size_t"

    @PYB11pure_virtual
    def Zero(self):
        "Set all element values equal to zero"
        return "void"
