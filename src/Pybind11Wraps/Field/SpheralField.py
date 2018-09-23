"""
Spheral Field module.

Provides the Field classes.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()
import copy

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Field/FieldBase.hh"',
            '"Field/Field.hh"',
            '"Field/FieldList.hh"',
            '"Field/FieldListSet.hh"',
            '"Utilities/FieldDataTypeTraits.hh"',
            '<vector>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# preamble
#-------------------------------------------------------------------------------
# preamble = ""
# for ndim in dims:
#     preamble += "typedef std::pair<NodeList<Dim<%(ndim)i>>*, std::string> pair_NodeList%(ndim)idptr_string;\n" % {"ndim": ndim}

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
FieldStorageType = PYB11enum(("ReferenceFields", "CopyFields"), export_values=True,
                             doc="Storage types for Fields in FieldLists.")

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from FieldBase import *
from Field import *

for ndim in (1,): #dims:

    #...........................................................................
    # FieldBase
    exec('''
FieldBase%(ndim)id = PYB11TemplateClass(FieldBase, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})

    #...........................................................................
    # non-numeric type fields.
    for (value, label) in (("Dim<%i>::FacetedVolume" % ndim,       "FacetedVolume"), 
                           ("std::vector<double>",                 "VectorDouble"),
                           ("std::vector<Dim<%i>::Vector>" % ndim, "VectorVector"),
                           ("std::vector<Dim<%i>::Tensor>" % ndim, "VectorSymTensor"),
                           ("std::vector<Dim<%i>::Tensor>" % ndim, "VectorSymTensor")):
        exec('''
%(label)sField%(ndim)sd = PYB11TemplateClass(Field, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    #...........................................................................
    # arithmetic fields
    for (value, label) in (("Dim<%i>::Vector" % ndim,           "Vector"),
                           ("Dim<%i>::Tensor" % ndim,           "Tensor"),
                           ("Dim<%i>::SymTensor" % ndim,        "SymTensor"),
                           ("Dim<%i>::ThirdRankTensor" % ndim,  "ThirdRankTensor"),
                           ("Dim<%i>::FourthRankTensor" % ndim, "FourthRankTensor"),
                           ("Dim<%i>::FifthRankTensor" % ndim,  "FifthRankTensor")):
        exec('''
%(label)sField%(ndim)sd = PYB11TemplateClass(ArithmeticField, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    #...........................................................................
    # Fully numeric fields
    for (value, label) in (("int",      "Int"),
                           ("uint64_t", "ULL"),
                           ("double",   "Scalar")):
        exec('''
%(label)sField%(ndim)sd = PYB11TemplateClass(MinMaxField, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})
