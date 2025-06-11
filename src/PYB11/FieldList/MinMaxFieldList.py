from PYB11Generator import *
import FieldList
import FieldListBase
import MinMaxFieldSpanList
from ArithmeticFieldList import *

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldList")
class MinMaxFieldList(FieldListBase.FieldListBase,
                      MinMaxFieldSpanList.MinMaxFieldSpanList):

    PYB11typedefs = """
    using FieldListType = FieldList<%(Dimension)s, %(Value)s>;
    using FieldType = Field<%(Dimension)s, %(Value)s>;
    using FieldSpanListType = FieldSpanList<%(Dimension)s, %(Value)s>;
    using FieldSpanType = FieldSpan<%(Dimension)s, %(Value)s>;
    using NodeListType = NodeList<%(Dimension)s>;
    using Scalar = %(Dimension)s::Scalar;
    using Vector = %(Dimension)s::Vector;
    using SymTensor = %(Dimension)s::SymTensor;
    using ViewType = typename FieldListType::ViewType;
"""

    def applyScalarMin(self):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self):
        "Enforce a double ceiling on the values of the Field."
        return

#-------------------------------------------------------------------------------
# Inject FieldList
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticFieldList, MinMaxFieldList)
