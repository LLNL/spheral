"""
Spheral FieldList module.

Provides the FieldList classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from FieldListBase import *
from FieldList import *
from FieldListSet import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Geometry/Dimension.hh"',
                  '"Field/FieldBase.hh"',
                  '"Field/Field.hh"',
                  '"Field/FieldList.hh"',
                  '"Field/FieldListSet.hh"',
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

    #...........................................................................
    # FieldListBase, FieldListSet
    exec('''
FieldListBase%(ndim)id = PYB11TemplateClass(FieldListBase, template_parameters="Dim<%(ndim)i>")
FieldListSet%(ndim)sd = PYB11TemplateClass(FieldListSet, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})

    #...........................................................................
    # non-numeric types
    for (value, label) in (("Dim<%i>::FacetedVolume" % ndim,       "FacetedVolume"), 
                           ("std::vector<int>",                    "VectorInt"),
                           ("std::vector<double>",                 "VectorDouble"),
                           ("std::vector<Dim<%i>::Vector>" % ndim, "VectorVector"),
                           ("std::vector<Dim<%i>::Tensor>" % ndim, "VectorSymTensor"),
                           ("std::vector<Dim<%i>::Tensor>" % ndim, "VectorSymTensor"),
                           ("std::vector<CellFaceFlag>",           "vector_of_CellFaceFlag"),
                           ("DomainNode<Dim<%i>>" % ndim,          "DomainNode"),
                           ("RKCoefficients<Dim<%i>>" % ndim,      "RKCoefficients")):
        exec('''
%(label)sFieldList%(ndim)sd = PYB11TemplateClass(FieldList, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    #...........................................................................
    # STL collections of Field types
    for value, label in (("int", "Int"),
                         ("double", "Scalar"),
                         ("Dim<%i>::Vector" % ndim, "Vector"),
                         ("Dim<%i>::Tensor" % ndim, "Tensor"),
                         ("Dim<%i>::SymTensor" % ndim, "SymTensor")):
        exec('''
vector_of_%(label)sFieldList%(ndim)id = PYB11_bind_vector("FieldList<%(Dimension)s, %(value)s>", opaque=True, local=False)
vector_of_%(label)sFieldListPtr%(ndim)id = PYB11_bind_vector("FieldList<%(Dimension)s, %(value)s>*", opaque=True, local=False)
vector_of_vector_of_%(label)sFieldList%(ndim)id = PYB11_bind_vector("std::vector<FieldList<%(Dimension)s, %(value)s>>", opaque=True, local=False)
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label,
       "Dimension" : "Dim<" + str(ndim) + ">"})
