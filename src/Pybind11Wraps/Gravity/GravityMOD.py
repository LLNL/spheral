"""
Spheral Gravity module.

Provides the implementations of gravitational physics packages in Spheral.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

from NBodyGravity import *
from TreeGravity import *

from GenericBodyForce import GenericBodyForce
PYB11import(GenericBodyForce, "SpheralPhysics")

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Physics/GenericBodyForce.hh"',
            '"Gravity/NBodyGravity.hh"',
            '"Gravity/TreeGravity.hh"',
            '"FileIO/FileIO.hh"',
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
for ndim in dims:
    exec('''
NBodyGravity%(ndim)id = PYB11TemplateClass(NBodyGravity, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

if 2 in dims:
    QuadTreeGravity = PYB11TemplateClass(TreeGravity, template_parameters="Dim<2>")

if 3 in dims:
    OctTreeGravity = PYB11TemplateClass(TreeGravity, template_parameters="Dim<3>")
