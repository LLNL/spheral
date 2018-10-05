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
# Instantiate our dimensional types
#-------------------------------------------------------------------------------
from SpheralFunctor import *
from KeyTraits import *
from Timer import *

ScalarScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "double"))
ScalarPairScalarFunctor = PYB11TemplateClass(SpheralFunctor, template_parameters=("double", "std::pair<double,double>"))

@PYB11template("Vector")
@PYB11implementation("[](std::vector<%(Vector)s>& positions) { %(Vector)s xmin, xmax; boundingBox(positions, xmin, xmax); return py::make_tuple(xmin, xmax); }")
def boundingBoxVec(positions = "std::vector<%(Vector)s>&"):
    "Minimum (axis-aligned) bounding box for a collection of %(Vector)s"
    return "py::tuple"

@PYB11template("Dimension")
@PYB11implementation("[](const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>& positions, const bool useGhosts) { %(Dimension)s::Vector xmin, xmax; boundingBox(positions, xmin, xmax, useGhosts); return py::make_tuple(xmin, xmax); }")
def boundingBoxFL(positions = "const FieldList<%(Dimension)s, typename %(Dimension)s::Vector>&",
                  useGhosts = "const bool"):
    "Minimum (axis-aligned) bounding box for a FieldList<%(Dimension)s::Vector>"
    return "py::tuple"

for ndim in dims:
    exec('''
# Functors
VectorScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "double"))
VectorVectorFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "%(Dimension)s::Vector"))
VectorPairScalarFunctor%(ndim)id = PYB11TemplateClass(SpheralFunctor, template_parameters=("%(Dimension)s::Vector", "std::pair<double,double>"))

# boundingVolumes
boundingBoxVec%(ndim)id = PYB11TemplateFunction(boundingBoxVec, template_parameters="%(Dimension)s::Vector", pyname="boundingBox")
boundingBoxFL%(ndim)id = PYB11TemplateFunction(boundingBoxFL, template_parameters="%(Dimension)s", pyname="boundingBox")
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

#-------------------------------------------------------------------------------
# packElement/unpackElement
#-------------------------------------------------------------------------------
@PYB11template("T")
def packElement(x = "const %(T)s&",
                buffer = "std::vector<char>&"):
    "Serialize a %(T)s onto a buffer of vector<char>"
    return "void"

@PYB11template("T")
@PYB11cppname("packElement")
def packElementVector(x = "const std::vector<%(T)s>&",
                      buffer = "std::vector<char>&"):
    "Serialize a vector<%(T)s> onto a buffer of vector<char>"
    return "void"

@PYB11template("T")
@PYB11implementation("[](const %(T)s& x) { std::vector<char> buf; packElement(x, buf); return std::string(buf.begin(), buf.end()); }")
def toString(x = "const %(T)s&"):
    return "std:string"

@PYB11template("T")
@PYB11implementation("[](const std::string& x) { %(T)s result; const std::vector<char> buf(x.begin(), x.end()); auto itr = buf.begin(); unpackElement(result, itr, buf.end()); return result; }")
def fromString(x = "const std::string&"):
    return "%(T)s"

for type, suffix in (("int", "Int"),
                     ("double", "Double"),
                     ("uint32_t", "Uint32_t"),
                     ("uint64_t", "Uint64_t"),
                     ("Dim<1>::FacetedVolume", "Box1d"),
                     ("Dim<2>::FacetedVolume", "Polygon"),
                     ("Dim<3>::FacetedVolume", "Polyhedron")):
    exec('''
packElement%(suffix)s       = PYB11TemplateFunction(packElement,       template_parameters="%(type)s",              pyname="packElement")
packElementVector%(suffix)s = PYB11TemplateFunction(packElementVector, template_parameters="%(type)s",              pyname="packElement")
toString%(suffix)s          = PYB11TemplateFunction(toString,          template_parameters="%(type)s",              pyname="toString")
toStringVector%(suffix)s    = PYB11TemplateFunction(toString,          template_parameters="std::vector<%(type)s>", pyname="toString")
to%(suffix)s                = PYB11TemplateFunction(fromString,        template_parameters="%(type)s")
toVector%(suffix)s          = PYB11TemplateFunction(fromString,        template_parameters="std::vector<%(type)s>")
''' % {"type"   : type,
       "suffix" : suffix})

packElementVecVecU     = PYB11TemplateFunction(packElementVector, template_parameters="std::vector<unsigned>",              pyname="packElement")
toStringVecVecU        = PYB11TemplateFunction(toString,          template_parameters="std::vector<std::vector<unsigned>>", pyname="toString")
toVectorVectorUnsigned = PYB11TemplateFunction(fromString,        template_parameters="std::vector<std::vector<unsigned>>")

#-------------------------------------------------------------------------------
# Geometry operations only available in 2D/3D.
#-------------------------------------------------------------------------------
@PYB11template("Vector")
def closestPointOnSegment(p = "const %(Vector)s&",
                          a0 = "const %(Vector)s&",
                          a1 = "const %(Vector)s&"):
    "Find the point on a line segment (a0,a1) closest to point (p)."
    return "%(Vector)s"

@PYB11template("Vector")
def pointPlaneDistance(point = "const %(Vector)s&",
                       origin = "const %(Vector)s&",
                       unitNormal = "const %(Vector)s&"):
    "Compute the distance between a point and a plane"
    return "double"

@PYB11template("Dimension")
def overlayRemapFields(boundaries = "const std::vector<Boundary<%(Dimension)s>*>&",
                       scalarDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::Scalar>*>&",
                       vectorDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::Vector>*>&",
                       tensorDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::Tensor>*>&",
                       symTensorDonorFields = "const std::vector<Field<%(Dimension)s, typename %(Dimension)s::SymTensor>*>&",
                       scalarAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::Scalar>*>&",
                       vectorAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::Vector>*>&",
                       tensorAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::Tensor>*>&",
                       symTensorAcceptorFields = "std::vector<Field<%(Dimension)s, typename %(Dimension)s::SymTensor>*>&"):
    """Use geometric clipping to remap a set of conserved fields.
Currently only works single NodeList -> single NodeList, no boundaries."""
    return "void"

for ndim in (x for x in dims if x in (2, 3)):
    exec('''
closestPointOnSegment%(ndim)id = PYB11TemplateFunction(closestPointOnSegment, template_parameters="%(Vector)s", pyname="closestPointOnSegment")
pointPlaneDistance%(ndim)id = PYB11TemplateFunction(pointPlaneDistance, template_parameters="%(Vector)s", pyname="pointPlaneDistance")
overlayRemapFields%(ndim)id = PYB11TemplateFunction(overlayRemapFields, template_parameters="Dim<%(ndim)i>", pyname="overlayRemapFields")
''' % {"ndim"   : ndim,
       "Vector" : "Dim<" + str(ndim) + ">::Vector"})
