from PYB11Generator import *
import FieldList
import FieldListBase
from ArithmeticFieldList import *

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldList")
class MinMaxFieldList(FieldListBase.FieldListBase):

    PYB11typedefs = """
    typedef FieldList<%(Dimension)s, %(Value)s> FieldListType;
    typedef Field<%(Dimension)s, %(Value)s> FieldType;
    typedef NodeList<%(Dimension)s> NodeListType;
    typedef %(Dimension)s::Scalar Scalar;
    typedef %(Dimension)s::Vector Vector;
    typedef %(Dimension)s::SymTensor SymTensor;
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
