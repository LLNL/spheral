"""
Spheral FieldSpanList module.

Provides the FieldSpanList classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#from FieldSpanListBase import *
from FieldSpanList import *

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
    FacetedVolume = f"{Dimension}::FacetedVolume"

#     #...........................................................................
#     # FieldSpanListBase
#     exec(f'''
# FieldSpanListBase{ndim}d = PYB11TemplateClass(FieldSpanListBase, template_parameters="{Dimension}")
# ''')

    #...........................................................................
    # FieldSpan -- non-numeric types 
    for (value, label) in (( FacetedVolume,                 "FacetedVolume"), 
                           ( "std::vector<int>",            "VectorInt"),
                           ( "std::vector<double>",         "VectorDouble"),
                           (f"std::vector<{Vector}>",       "VectorVector"),
                           (f"std::vector<{Tensor}>",       "VectorTensor"),
                           (f"std::vector<{SymTensor}>",    "VectorSymTensor"),
                           ( "std::vector<CellFaceFlag>",   "vector_of_CellFaceFlag"),
                           (f"DomainNode<{Dimension}>",     "DomainNode"),
                           (f"RKCoefficients<{Dimension}>", "RKCoefficients")):
        exec(f'''
{label}FieldSpanList{ndim}d = PYB11TemplateClass(FieldSpanList, template_parameters=("{Dimension}", "{value}"))
''')
