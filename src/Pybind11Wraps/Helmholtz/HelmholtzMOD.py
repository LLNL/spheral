"""
Spheral Helmholtz module.

Provides the Helmholtz equation of state
"""

from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from HelmholtzEquationOfState import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"Material/PhysicalConstants.hh"',
                  '"Material/EquationOfState.hh"',
                  '"Material/HelmholtzEquationOfState.hh"']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
PYB11namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our dimensional types
#-------------------------------------------------------------------------------
for ndim in dims:
    exec('''
HelmholtzEquationOfState%(ndim)id = PYB11TemplateClass(HelmholtzEquationOfState, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
