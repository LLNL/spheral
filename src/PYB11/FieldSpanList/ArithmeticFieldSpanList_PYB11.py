"""
Spheral ArithmeticFieldSpanList module.

Provides the ArithmeticFieldSpanList classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from ArithmeticFieldSpanList import *
from MinMaxFieldSpanList import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldSpanBase.hh"',
                  '"Field/FieldSpan.hh"',
                  '"Field/FieldSpanList.hh"',
                  '"Utilities/FieldDataTypeTraits.hh"',
                  '"Utilities/DomainNode.hh"',
                  '"Geometry/CellFaceFlag.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    Vector = f"{Dimension}::Vector"
    Tensor = f"{Dimension}::Tensor"
    SymTensor = f"{Dimension}::SymTensor"

    #...........................................................................
    # arithmetic fields
    for (value, label) in (("int",            "Int"),
                           ("unsigned",       "Unsigned"),
                           ("uint64_t",       "ULL"),
                           (Vector,           "Vector"),
                           (Tensor,           "Tensor")):
        exec(f'''
{label}FieldSpanList{ndim}d = PYB11TemplateClass(ArithmeticFieldSpanList, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # A few fields can apply the min/max with a scalar additionally
    for (value, label) in (("double",   "Scalar"),
                           (SymTensor,  "SymTensor")):
        exec(f'''
{label}FieldSpanList{ndim}d = PYB11TemplateClass(MinMaxFieldSpanList, template_parameters=("{Dimension}", "{value}"))
''')
