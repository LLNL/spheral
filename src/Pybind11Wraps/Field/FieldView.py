from PYB11Generator import *
from Field import Field

#-------------------------------------------------------------------------------
# FieldView
#-------------------------------------------------------------------------------
@PYB11template("Dimension", "Value")
@PYB11module("SpheralField")
class FieldView():
    "FieldView class"

    PYB11typedefs = """
    typedef FieldView<%(Dimension)s, %(Value)s> FieldViewType;
    typedef Field<%(Dimension)s, %(Value)s> FieldType;
"""

