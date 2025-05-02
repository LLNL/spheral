import inspect
from PYB11Generator import *
from FieldSpanBase import *
from ArithmeticFieldSpan import *

#-------------------------------------------------------------------------------
# Add min/max operations to a FieldSpan
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11pycppname("FieldSpan")
class MinMaxFieldSpan(FieldSpanBase):

    PYB11typedefs = """
  using SelfType = FieldSpan<%(Dimension)s, %(Value)s>;
  using Scalar = typename SelfType::Scalar;
  using ScalarFieldSpan = FieldSpan<%(Dimension)s, Scalar>;
"""

    def applyScalarMin(self):
        "Enforce a double floor on the values of the FieldSpan."
        return

    def applyScalarMax(self):
        "Enforce a double ceiling on the values of the FieldSpan."
        return

#-------------------------------------------------------------------------------
# Inject base field methods
#-------------------------------------------------------------------------------
PYB11inject(ArithmeticFieldSpan, MinMaxFieldSpan)
