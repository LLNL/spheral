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
            '"Field/Field.hh"',
            '"Field/FieldList.hh"',
            '"FileIO/FileIO.hh"',
            '<vector>',
            '<string>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------
@PYB11template("Boundary1", "Boundary2")
@PYB11returnpolicy("reference")
@PYB11implementation("[](%(Boundary1)s* boundPtr) { return dynamic_cast<%(Boundary2)s*>(boundPtr); }")
def dynamicCastBoundary(boundPtr = "%(Boundary1)s*"):
    "Dynamic cast %(Boundary1)s* -> %(Boundary2)s*"
    return "%(Boundary2)s*"

#-------------------------------------------------------------------------------
# Do our dimension dependent instantiations.
#-------------------------------------------------------------------------------
from Boundary import *
from PlanarBoundary import *
from ReflectingBoundary import *
from RigidBoundary import *
from PeriodicBoundary import *
from ConstantVelocityBoundary import *
from ConstantXVelocityBoundary import *
from ConstantYVelocityBoundary import *
from ConstantZVelocityBoundary import *
from ConstantRVelocityBoundary import *
from ConstantBoundary import *
from CRKSPHVoidBoundary import *

for ndim in dims:
    exec('''
Boundary%(ndim)id = PYB11TemplateClass(Boundary, template_parameters="%(Dimension)s")
PlanarBoundary%(ndim)id = PYB11TemplateClass(PlanarBoundary, template_parameters="%(Dimension)s")
ReflectingBoundary%(ndim)id = PYB11TemplateClass(ReflectingBoundary, template_parameters="%(Dimension)s")
RigidBoundary%(ndim)id = PYB11TemplateClass(RigidBoundary, template_parameters="%(Dimension)s")
PeriodicBoundary%(ndim)id = PYB11TemplateClass(PeriodicBoundary, template_parameters="%(Dimension)s")
ConstantVelocityBoundary%(ndim)id = PYB11TemplateClass(ConstantVelocityBoundary, template_parameters="%(Dimension)s")
ConstantXVelocityBoundary%(ndim)id = PYB11TemplateClass(ConstantXVelocityBoundary, template_parameters="%(Dimension)s")
ConstantBoundary%(ndim)id = PYB11TemplateClass(ConstantBoundary, template_parameters="%(Dimension)s")
CRKSPHVoidBoundary%(ndim)id = PYB11TemplateClass(CRKSPHVoidBoundary, template_parameters="%(Dimension)s")

vector_of_Boundary%(ndim)id = PYB11_bind_vector("Boundary<%(Dimension)s>*", opaque=True)

dynamicCastBoundary%(ndim)id = PYB11TemplateFunction(dynamicCastBoundary, template_parameters = ("Boundary<%(Dimension)s>", "PlanarBoundary<%(Dimension)s>"))
''' % {"ndim"      : ndim,
       "Dimension" : ("Dim<" + str(ndim) +">")})

if 2 in dims:
    ConstantYVelocityBoundary2d = PYB11TemplateClass(ConstantYVelocityBoundary, template_parameters="Dim<2>")
    from AxisBoundaryRZ import *

if 3 in dims:
    ConstantYVelocityBoundary2d = PYB11TemplateClass(ConstantYVelocityBoundary, template_parameters="Dim<2>")
    ConstantYVelocityBoundary3d = PYB11TemplateClass(ConstantYVelocityBoundary, template_parameters="Dim<3>")
    ConstantZVelocityBoundary3d = PYB11TemplateClass(ConstantZVelocityBoundary, template_parameters="Dim<3>")
    from CylindricalBoundary import *
    from SphericalBoundary import *

