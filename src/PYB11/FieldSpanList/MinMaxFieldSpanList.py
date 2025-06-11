from PYB11Generator import *
import FieldSpanList
from ArithmeticFieldSpanList import *

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldSpanList")
class MinMaxFieldSpanList:

    PYB11typedefs = """
    using FieldType = Field<%(Dimension)s, %(Value)s>;
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldSpanType = FieldSpan<%(Dimension)s, %(Value)s>;
    using FieldSpanListType = FieldSpanList<%(Dimension)s, %(Value)s>;
    using Scalar = %(Dimension)s::Scalar;
    using Vector = %(Dimension)s::Vector;
    using SymTensor = %(Dimension)s::SymTensor;
"""

    def applyMin(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s floor on the values of the Field."
        return

    def applyMax(self, rhs="const %(Value)s&"):
        "Enforce a %(Value)s ceiling on the values of the Field."
        return

    def applyScalarMin(self, rhs="const Scalar"):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self, rhs="const Scalar"):
        "Enforce a double ceiling on the values of the Field."
        return

#-------------------------------------------------------------------------------
# Inject FieldSpanList
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticFieldSpanList, MinMaxFieldSpanList)
