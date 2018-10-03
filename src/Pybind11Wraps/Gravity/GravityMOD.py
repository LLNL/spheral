"""
Spheral Gravity module.

Provides the implementations of gravitational physics packages in Spheral.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Physics/GenericBodyForce.hh"',
            '"Gravity/NBodyGravity.hh"',
            '"Gravity/TreeGravity.hh"',
            '<vector>',
            '<string>',
            '<iterator>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Enums
#-------------------------------------------------------------------------------
GravityTimeStepType = PYB11enum(("AccelerationRatio", "DynamicalTime"), export_values=True,
                                doc="Algorithm choice for choosing timesteps for TreeGravity")

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from NBodyGravity import *
from TreeGravity import *

for ndim in dims:
    exec('''
NBodyGravity%(ndim)id = PYB11TemplateClass(NBodyGravity, template_parameters="%(Dimension)s")
TreeGravity%(ndim)id = PYB11TemplateClass(TreeGravity, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
