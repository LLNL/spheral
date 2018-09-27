"""
Spheral Boundary module.

Provides the Boundary base class and many boundary implementations.
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Boundary/Boundary.hh"',
            '"Boundary/PlanarBoundary.hh"',
            '"Boundary/ReflectingBoundary.hh"',
            '"Boundary/RigidBoundary.hh"',
            '"Boundary/PeriodicBoundary.hh"',
            '"Boundary/ConstantVelocityBoundary.hh"',
            '"Boundary/ConstantXVelocityBoundary.hh"',
            '"Boundary/ConstantYVelocityBoundary.hh"',
            '"Boundary/ConstantZVelocityBoundary.hh"',
            '"Boundary/ConstantRVelocityBoundary.hh"',
            '"Boundary/ConstantBoundary.hh"',
            '"Boundary/SphericalBoundary.hh"',
            '"Boundary/CylindricalBoundary.hh"',
            '"Boundary/AxialSymmetryBoundary.hh"',
            '"Boundary/AxisBoundaryRZ.hh"',
            '"Boundary/CRKSPHVoidBoundary.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from Boundary import *
from PlanarBoundary import *

for ndim in dims:
    exec('''
Boundary%(ndim)id = PYB11TemplateClass(Boundary, template_parameters="Dim<%(ndim)i>")
PlanarBoundary%(ndim)id = PYB11TemplateClass(PlanarBoundary, template_parameters="Dim<%(ndim)i>")
''' % {"ndim" : ndim})


