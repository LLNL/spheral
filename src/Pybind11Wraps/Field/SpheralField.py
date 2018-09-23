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

# ArithmeticField = copy.deepcopy(Field)
# addFieldArithmeticOperations(ArithmeticField)
# ArithmeticField.PYB11cppname = ArithmeticField.PYB11pyname = "Field"

# NumericField = copy.deepcopy(Field)
# addFieldArithmeticOperations(NumericField)
# addFieldMinMaxOperations(Field)
# NumericField.PYB11cppname = NumericField.PYB11pyname = "Field"

for ndim in (1,): #dims:
    exec('''
FieldBase%(ndim)id = PYB11TemplateClass(FieldBase, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})

    # First the non-numeric type fields.
    for (value, label) in (("std::vector<double>", "VectorDouble"),):
        exec('''
%(label)sField%(ndim)sd = PYB11TemplateClass(Field, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

# Arithmetic only fields
for ndim in (1,): #dims:
    for (value, label) in (("Dim<%i>::Vector" % ndim, "Vector"),):
        exec('''
%(label)sField%(ndim)sd = PYB11TemplateClass(ArithmeticField, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

# Fully numeric fields
for ndim in (1,): #dims:
    for (value, label) in (("double", "Scalar"),
                           ("int", "Int")):
        exec('''
%(label)sField%(ndim)sd = PYB11TemplateClass(MinMaxField, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})
