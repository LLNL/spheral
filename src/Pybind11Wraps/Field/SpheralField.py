"""
Spheral Field module.

Provides the Field classes.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

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
    exec('''
FieldBase%(ndim)id = PYB11TemplateClass(FieldBase, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})

    # Arithmetic fields
    for (value, label) in (("double", "Scalar"),
                           ("int", "Int")):
        exec('''
addFieldArithmeticOperations(Field)
addFieldMinMaxOperations(Field)
%(label)sField%(ndim)sd = PYB11TemplateClass(Field, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    # Arithmetic only fields
    for (value, label) in (("Dim<%i>::Vector" % ndim, "Vector"),):
        exec('''
Field = None
from Field import Field
addFieldArithmeticOperations(Field)
%(label)sField%(ndim)sd = PYB11TemplateClass(Field, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    # Non-numeric fields
    for (value, label) in (("std::vector<double>", "VectorDouble"),):
        exec('''
Field = None
from Field import Field
%(label)sField%(ndim)sd = PYB11TemplateClass(Field, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})
