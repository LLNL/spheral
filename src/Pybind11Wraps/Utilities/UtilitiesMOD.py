"""
Spheral Utilities module.

A bunch of utility methods used throughout Spheral.  Unfortunately this has become
a bit of a grab bag of math, geometry, infrastructure, and assorted functions with
no real relation.  Should probably revisit and categorize this stuff to other
modules more effectively.a
"""

from PYB11Generator import *
from spheralDimensions import *
dims = spheralDimensions()

#-------------------------------------------------------------------------------
# Includes
#-------------------------------------------------------------------------------
includes = ['"Geometry/Dimension.hh"',
            '"Utilities/packElement.hh"',
            '"boost/math/special_functions/legendre.hpp"',
            '"Utilities/Functors.hh"',
            '"Utilities/erff.hh"',
            '"Utilities/newtonRaphson.hh"',
            '"Utilities/simpsonsIntegration.hh"',
            '"Utilities/globalNodeIDs.hh"',
            '"Utilities/rotationMatrix.hh"',
            '"Utilities/iterateIdealH.hh"',
            '"Utilities/mortonOrderIndices.hh"',
            '"Utilities/peanoHilbertOrderIndices.hh"',
            '"Utilities/boundingBox.hh"',
            '"Utilities/globalBoundingVolumes.hh"',
            '"Utilities/testBoxIntersection.hh"',
            '"Utilities/lineSegmentIntersections.hh"',
            '"Utilities/pointDistances.hh"',
            '"Utilities/segmentIntersectEdges.hh"',
            '"Utilities/integrateThroughMeshAlongSegment.hh"',
            '"Utilities/numberDensity.hh"',
            '"Utilities/planarReflectingOperator.hh"',
            '"Utilities/KeyTraits.hh"',
            '"Utilities/pointInPolygon.hh"',
            '"Utilities/pointOnPolygon.hh"',
            '"Utilities/pointInPolyhedron.hh"',
            '"Utilities/pointOnPolyhedron.hh"',
            '"Utilities/refinePolyhedron.hh"',
            '"Utilities/overlayRemapFields.hh"',
            '"Utilities/computeShepardsInterpolation.hh"',
            '"Utilities/Timer.hh"',
            '<vector>',
            '<string>',
            '<iterator>']

#-------------------------------------------------------------------------------
# Namespaces
#-------------------------------------------------------------------------------
namespaces = ["Spheral"]

#-------------------------------------------------------------------------------
# Instantiate our types
#-------------------------------------------------------------------------------
from SpheralFunctor import *

ScalarScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "double"))
ScalarPairScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "std::pair<double,double>"))

for ndim in dims:
    exec('''
VectorScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "double"))
VectorVectorFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "%(Dimension)s::Vector"))
VectorPairScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "std::pair<double,double>"))
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})
