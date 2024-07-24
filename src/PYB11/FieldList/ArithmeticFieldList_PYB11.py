"""
Spheral ArithmeticFieldList module.

Provides the ArithmeticFieldList classes.
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from ArithmeticFieldList import *
from MinMaxFieldList import *

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
    # arithmetic fields
    for (value, label) in (("int",                              "Int"),
                           ("unsigned",                         "Unsigned"),
                           ("uint64_t",                         "ULL"),
                           ("Dim<%i>::Vector" % ndim,           "Vector"),
                           ("Dim<%i>::Tensor" % ndim,           "Tensor"),
                           ("Dim<%i>::ThirdRankTensor" % ndim,  "ThirdRankTensor"),
                           ("Dim<%i>::FourthRankTensor" % ndim, "FourthRankTensor"),
                           ("Dim<%i>::FifthRankTensor" % ndim,  "FifthRankTensor")):
        exec('''
%(label)sFieldList%(ndim)sd = PYB11TemplateClass(ArithmeticFieldList, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})

    #...........................................................................
    # A few fields can apply the min/max with a scalar additionally
    for (value, label) in (("double",                     "Scalar"),
                           ("Dim<%i>::SymTensor" % ndim,  "SymTensor")):
        exec('''
%(label)sFieldList%(ndim)sd = PYB11TemplateClass(MinMaxFieldList, template_parameters=("Dim<%(ndim)i>", "%(value)s"))
''' % {"ndim" : ndim,
       "value" : value,
       "label" : label})
