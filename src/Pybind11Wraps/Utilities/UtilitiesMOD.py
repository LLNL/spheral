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
from KeyTraits import *
from Timer import *

ScalarScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "double"))
ScalarPairScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "std::pair<double,double>"))

for ndim in dims:
    exec('''
VectorScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "double"))
VectorVectorFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "%(Dimension)s::Vector"))
VectorPairScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "std::pair<double,double>"))
''' % {"ndim"      : ndim,
       "Dimension" : "Dim<" + str(ndim) + ">"})

#-------------------------------------------------------------------------------
# Module functions
#-------------------------------------------------------------------------------
def erff(x = "double"):
    "You know, the error function."
    return "double"

@PYB11cppname("newtonRaphson<const PythonBoundFunctors::SpheralFunctor<double, std::pair<double, double>>>")
def newtonRaphsonFindRoot(function = "const PythonBoundFunctors::SpheralFunctor<double, std::pair<double, double>>&",
                          x1 = "double",
                          x2 = "double",
                          xaccuracy = ("double", "1.0e-15"),
                          yaccuracy = ("double", "1.0e-15"),
                          maxIterations = ("int", "100")):
    """Newton-Raphson root finder.
Finds a root of 'function' in the range (x1, x2)"""
    return "double"

@PYB11cppname("simpsonsIntegration<const PythonBoundFunctors::SpheralFunctor<double, double>, double, double>")
def simpsonsIntegrationDouble(function = "const PythonBoundFunctors::SpheralFunctor<double, double>&",
                              x0 = "double",
                              x1 = "double",
                              numBins = "unsigned"):
    "Numerically integrate 'function' in the range (x0, x1) via Simpsons rule"
    return "double"
