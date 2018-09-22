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
from FieldBase import FieldBase

for ndim in (1,): #dims:
    exec('''
FieldBase%(ndim)id = PYB11TemplateClass(FieldBase, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})
