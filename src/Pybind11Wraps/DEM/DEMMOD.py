"""
Spheral DEM module.

Provides implementations of DEM
"""
from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from DEMBase import *
from ContactModelBase import *
#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"DEM/DEMBase.hh"',
                  '"DEM/ContactModelBase.hh"',
                  '"FileIO/FileIO.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]


#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
DEMBase%(ndim)id = PYB11TemplateClass(DEMBase, template_parameters="%(Dimension)s")
ContactModelBase%(ndim)id = PYB11TemplateClass(ContactModelBase, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

