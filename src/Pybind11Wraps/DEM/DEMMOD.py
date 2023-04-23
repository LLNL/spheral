"""
Spheral DEM module.

Provides implementations of DEM
"""
from PYB11Generator import *
from SpheralCommon import *
from spheralDimensions import *
dims = spheralDimensions()

from DEMBase import *
from LinearSpringDEM import *
from SolidBoundaries import *

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
PYB11includes += ['"DEM/DEMBase.hh"',
                  '"DEM/DEMFieldNames.hh"',
                  '"DEM/LinearSpringDEM.hh"',
                  '"DEM/SolidBoundary/SolidBoundary.hh"',
                  '"DEM/SolidBoundary/PlanarWall.hh"',
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
LinearSpringDEM%(ndim)id = PYB11TemplateClass(LinearSpringDEM, template_parameters="%(Dimension)s")
SolidBoundary%(ndim)id = PYB11TemplateClass(SolidBoundary, template_parameters="%(Dimension)s")
PlanarWall%(ndim)id = PYB11TemplateClass(PlanarWall, template_parameters="%(Dimension)s")
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

#-------------------------------------------------------------------------------
# expose our field names
#-------------------------------------------------------------------------------
class DEMFieldNames:
    particleRadius = PYB11readonly(static=True, returnpolicy="copy")
    compositeParticleIndex = PYB11readonly(static=True, returnpolicy="copy")
    angularVelocity = PYB11readonly(static=True, returnpolicy="copy")
    uniqueIndices = PYB11readonly(static=True, returnpolicy="copy")
    isActiveContact = PYB11readonly(static=True, returnpolicy="copy")
    neighborIndices = PYB11readonly(static=True, returnpolicy="copy")
    shearDisplacement = PYB11readonly(static=True, returnpolicy="copy")
    rollingDisplacement = PYB11readonly(static=True, returnpolicy="copy")
    torsionalDisplacement = PYB11readonly(static=True, returnpolicy="copy")
    equilibriumOverlap = PYB11readonly(static=True, returnpolicy="copy")
    solidBoundaries = PYB11readonly(static=True, returnpolicy="copy")