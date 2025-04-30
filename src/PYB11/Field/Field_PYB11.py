"""
Spheral Field module.

Provides the Field classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from FieldBase import *
from Field import *
from ArithmeticField import *
from MinMaxField import *
from FieldSpanBase import *
from FieldSpan import *
from ArithmeticFieldSpan import *
from MinMaxFieldSpan import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldBase.hh"',
                  '"Field/Field.hh"',
                  '"Field/FieldList.hh"',
                  '"Field/FieldListSet.hh"',
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
# Enums
#-------------------------------------------------------------------------------
FieldStorageType = PYB11enum(("ReferenceFields", "CopyFields"), export_values=True,
                             doc="Storage types for Fields in FieldLists.")

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
    # FieldBase
    exec(f'''
FieldBase{ndim}d = PYB11TemplateClass(FieldBase, template_parameters="{Dimension}")
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
{label}Field{ndim}d = PYB11TemplateClass(Field, template_parameters=("{Dimension}", "{value}"))
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
{label}Field{ndim}d = PYB11TemplateClass(ArithmeticField, template_parameters=("{Dimension}", "{value}"))
{label}FieldSpan{ndim}d = PYB11TemplateClass(ArithmeticFieldSpan, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # A few fields can apply the min/max with a scalar addtionally
    for (value, label) in ((Scalar,     "Scalar"),
                           (SymTensor,  "SymTensor")):
        exec(f'''
{label}Field{ndim}d = PYB11TemplateClass(MinMaxField, template_parameters=("{Dimension}", "{value}"))
{label}FieldSpan{ndim}d = PYB11TemplateClass(MinMaxFieldSpan, template_parameters=("{Dimension}", "{value}"))
''')

    #...........................................................................
    # STL collections of Field types
    for value, label in (("int",     "Int"),
                         (Scalar,    "Scalar"),
                         (Vector,    "Vector"),
                         (Tensor,    "Tensor"),
                         (SymTensor, "SymTensor")):
        exec(f'''
vector_of_{label}Field{ndim}d = PYB11_bind_vector("Field<{Dimension}, {value}>", opaque=True, local=False)
vector_of_{label}FieldPtr{ndim}d = PYB11_bind_vector("Field<{Dimension}, {value}>*", opaque=True, local=False)
vector_of_vector_of_{label}Field{ndim}d = PYB11_bind_vector("std::vector<Field<{Dimension}, {value}>>", opaque=True, local=False)
''')
