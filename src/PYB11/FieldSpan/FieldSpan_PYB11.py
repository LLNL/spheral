"""
Spheral FieldSpan module.

Provides the FieldSpan classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from FieldSpanBase import *
from FieldSpan import *
from ArithmeticFieldSpan import *
from MinMaxFieldSpan import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldSpanBase.hh"',
                  '"Field/FieldSpan.hh"',
                  '"Utilities/FieldDataTypeTraits.hh"',
                  '"Utilities/DomainNode.hh"',
                  '"Geometry/CellFaceFlag.hh"',
                  '<vector>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
for ndim in dims:
    Dimension = f"Dim<{ndim}>"
    Scalar = f"{Dimension}::Scalar"
    Vector = f"{Dimension}::Vector"
    Tensor = f"{Dimension}::Tensor"
    SymTensor = f"{Dimension}::SymTensor"
    ThirdRankTensor = f"{Dimension}::ThirdRankTensor"
    FourthRankTensor = f"{Dimension}::FourthRankTensor"
    FifthRankTensor = f"{Dimension}::FifthRankTensor"

    #...........................................................................
    # FieldSpanBase
    exec(f'''
FieldSpanBase{ndim}d = PYB11TemplateClass(FieldSpanBase, template_parameters="{Dimension}")
''')

    #...........................................................................
    # non-numeric types
    for (value, label) in ((f"{Dimension}::FacetedVolume",       "FacetedVolume"),
                           ( "std::vector<int>",                 "VectorInt"),
                           ( "std::vector<double>",              "VectorDouble"),
                           (f"std::vector<{Vector}>",            "VectorVector"),
                           (f"std::vector<{Tensor}>",            "VectorSymTensor"),
                           (f"std::vector<{Tensor}>",            "VectorSymTensor"),
                           ( "std::vector<CellFaceFlag>",        "vector_of_CellFaceFlag"),
                           (f"DomainNode<{Dimension}>",          "DomainNode"),
                           (f"RKCoefficients<{Dimension}>",      "RKCoefficients")):
        exec(f'''
{label}FieldSpan{ndim}d = PYB11TemplateClass(FieldSpan, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # arithmetic fields
    for (value, label) in (("int",            "Int"),
                           ("unsigned",       "Unsigned"),
                           ("uint64_t",       "ULL"),
                           (Vector,           "Vector"),
                           (Tensor,           "Tensor"),
                           (ThirdRankTensor,  "ThirdRankTensor"),
                           (FourthRankTensor, "FourthRankTensor"),
                           (FifthRankTensor,  "FifthRankTensor")):
        exec(f'''
{label}FieldSpan{ndim}d = PYB11TemplateClass(ArithmeticFieldSpan, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # A few fields can apply the min/max with a scalar addtionally
    for (value, label) in ((Scalar,     "Scalar"),
                           (SymTensor,  "SymTensor")):
        exec(f'''
{label}FieldSpan{ndim}d = PYB11TemplateClass(MinMaxFieldSpan, template_parameters=("{Dimension}", "{value}"))
''')
