"""
Spheral Field STL module.

Provides the STL containers for Field classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldBase.hh"',
                  '"Field/Field.hh"',
                  '"Field/FieldList.hh"',
                  '"Field/FieldListSet.hh"',
                  '"Utilities/FieldDataTypeTraits.hh"',
                  '<vector>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
for ndim in dims:

    #...........................................................................
    # STL collections of Field types
    for value, label in (("int", "Int"),
                         ("double", "Scalar"),
                         ("Dim<%i>::Vector" % ndim, "Vector"),
                         ("Dim<%i>::Tensor" % ndim, "Tensor"),
                         ("Dim<%i>::SymTensor" % ndim, "SymTensor")):
        exec('''
vector_of_%(label)sField%(ndim)id = PYB11_bind_vector("Field<%(Dimension)s, %(value)s>", opaque=True, local=False)
vector_of_%(label)sFieldPtr%(ndim)id = PYB11_bind_vector("Field<%(Dimension)s, %(value)s>*", opaque=True, local=False)
vector_of_vector_of_%(label)sField%(ndim)id = PYB11_bind_vector("std::vector<Field<%(Dimension)s, %(value)s>>", opaque=True, local=False)

vector_of_%(label)sFieldList%(ndim)id = PYB11_bind_vector("FieldList<%(Dimension)s, %(value)s>", opaque=True, local=False)
vector_of_%(label)sFieldListPtr%(ndim)id = PYB11_bind_vector("FieldList<%(Dimension)s, %(value)s>*", opaque=True, local=False)
vector_of_vector_of_%(label)sFieldList%(ndim)id = PYB11_bind_vector("std::vector<FieldList<%(Dimension)s, %(value)s>>", opaque=True, local=False)
''' % {"ndim": ndim,
       "Dimension" : "Dim<%i>" % ndim,
       "value" : value,
       "label" : label})
