import inspect
from PYB11Generator import *
from FieldBase import *
from ArithmeticField import *

#-------------------------------------------------------------------------------
# Add min/max operations to a Field
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("Field")
class MinMaxField(FieldBase):

    PYB11typedefs = """
  typedef Field<%(Dimension)s, %(Value)s> FieldType;
"""

    def applyScalarMin(self):
        "Enforce a double floor on the values of the Field."
        return

    def applyScalarMax(self):
        "Enforce a double ceiling on the values of the Field."
        return


@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldView")
class MinMaxFieldView(FieldView):
    PYB11typedefs = """
    typedef FieldView<%(Dimension)s, %(Value)s> FieldViewType;
    typedef Field<%(Dimension)s, %(Value)s> FieldType;
"""

#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticField, MinMaxField)
PYB11inject(ArithmeticFieldView, MinMaxFieldView)
