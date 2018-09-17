"""
This module provides the fundamental Spheral Geometry types (Vector, Tensor, 
etc.) and associated methods such as products and eigenvalues.
"""

from PYB11Decorators import *
from PYB11STLmethods import *
from PYB11property import *
from PYB11class import *
from PYB11function import *
import types

# Define some useful type collections we're going to be wrapping in this module.
geomtypes = ["Vector", "Tensor", "SymTensor", "ThirdRankTensor", "FourthRankTensor", "FifthRankTensor", "FacetedVolume"]

# The code preamble
preamble = """
using namespace Spheral;

"""
for ndim in (1, 2, 3):
    preamble += "typedef GeomPlane<Dim<%i>> Plane%id;\n" % (ndim, ndim)
    for gtype in geomtypes:
        preamble += "typedef Dim<%i>::%s %s%id;\n" % (ndim, gtype, gtype, ndim)

# Include files
includes = ['"Geometry/Dimension.hh"',
            '"Geometry/GeomVector.hh"',
            '"Geometry/Geom3Vector.hh"',
            '"Geometry/GeomTensor.hh"',
            '"Geometry/GeomSymmetricTensor.hh"',
            '"Geometry/GeomThirdRankTensor.hh"',
            '"Geometry/GeomFourthRankTensor.hh"',
            '"Geometry/GeomFifthRankTensor.hh"',
            '"Geometry/EigenStruct.hh"',
            '"Geometry/computeEigenValues.hh"',
            '"Geometry/GeomPlane.hh"',
            '"Geometry/GeomPolygon.hh"',
            '"Geometry/GeomPolyhedron.hh"',
            '"Geometry/GeomFacet2d.hh"',
            '"Geometry/GeomFacet3d.hh"',
            '"Geometry/invertRankNTensor.hh"',
            '"Geometry/innerProduct.hh"',
            '"Geometry/outerProduct.hh"',
            '"Geometry/innerDoubleProduct.hh"',
            '"Geometry/aggregateFacetedVolumes.hh"',
            '"Field/Field.hh"',
            '"Utilities/DataTypeTraits.hh"',

            '<vector>',
            '<sstream>']

# STL containers
for gtype in geomtypes + ["Plane"]:
    for suffix in ("1d", "2d", "3d"):
        element = gtype + suffix
        exec('vector_of_%s = PYB11_bind_vector("%s", opaque=True)' % (element, element))
vector_of_Facet2d = PYB11_bind_vector("GeomFacet2d", opaque=True)
vector_of_Facet3d = PYB11_bind_vector("GeomFacet3d", opaque=True)

# Get the objects wrapped in other files.
from Vector import Vector1d, Vector2d, Vector3d
from Tensor import Tensor1d, Tensor2d, Tensor3d
from SymTensor import SymTensor1d, SymTensor2d, SymTensor3d
from ThirdRankTensor import ThirdRankTensor1d, ThirdRankTensor2d, ThirdRankTensor3d
from FourthRankTensor import FourthRankTensor1d, FourthRankTensor2d, FourthRankTensor3d
from FifthRankTensor import FifthRankTensor1d, FifthRankTensor2d, FifthRankTensor3d
from EigenStruct import EigenStruct1d, EigenStruct2d, EigenStruct3d
from Plane import Plane1d, Plane2d, Plane3d

#-------------------------------------------------------------------------------
# invertRankNTensor template
#-------------------------------------------------------------------------------
@PYB11template("TensorType")
def invertRankNTensor(tensor = "const %(TensorType)s&"):
    "Compute the inverse of a tensor."
    return "%(TensorType)s"

invertRankNTensor1 = PYB11TemplateFunction(invertRankNTensor,
                                           template_parameters = "Dim<1>::Tensor",
                                           pyname = "invertRankNTensor")
invertRankNTensor2 = PYB11TemplateFunction(invertRankNTensor,
                                           template_parameters = "Dim<1>::SymTensor",
                                           pyname = "invertRankNTensor")
invertRankNTensor3 = PYB11TemplateFunction(invertRankNTensor,
                                           template_parameters = "Dim<1>::FourthRankTensor",
                                           pyname = "invertRankNTensor")

#-------------------------------------------------------------------------------
# computeEigenValues
#-------------------------------------------------------------------------------
@PYB11template("Dim")
def computeEigenValues(field = "const Field<%(Dim)s, %(Dim)s::SymTensor>&",
                       eigenValues = "const Field<%(Dim)s, %(Dim)s::Vector>&",
                       eigenVectors = "const Field<%(Dim)s, %(Dim)s::Tensor>&"):
    "Compute the eigenvalues for a field of symmetric tensors."
    return "void"

computeEigenValues1 = PYB11TemplateFunction(computeEigenValues,
                                            template_parameters = "Dim<1>",
                                            pyname = "computeEigenValues")
computeEigenValues2 = PYB11TemplateFunction(computeEigenValues,
                                            template_parameters = "Dim<2>",
                                            pyname = "computeEigenValues")
computeEigenValues3 = PYB11TemplateFunction(computeEigenValues,
                                            template_parameters = "Dim<3>",
                                            pyname = "computeEigenValues")

#-------------------------------------------------------------------------------
# Inner product (with a double)
#-------------------------------------------------------------------------------
@PYB11template("ValueType")
def innerProductScalar(A = "const double&",
                       B = "const %(ValueType)s&"):
    "Inner product with a scalar."
    return "%(ValueType)s"

@PYB11template("ValueType")
def innerProductScalarR(A = "const %(ValueType)s&",
                        B = "const double&"):
    "Inner product with a scalar."
    return "%(ValueType)s"

i = 0
for VT in ("Vector", "Tensor", "SymTensor", "ThirdRankTensor", "FourthRankTensor", "FifthRankTensor"):
    for ndim in (1, 2, 3):
        exec("""
innerProductScalar%(i)i = PYB11TemplateFunction(innerProductScalar,
                                                template_parameters = "Dim<%(ndim)i>::%(VT)s",
                                                pyname = "innerProduct",
                                                cppname = "innerProduct<Dim<%(ndim)i>::%(VT)s>")
innerProductScalar%(j)i = PYB11TemplateFunction(innerProductScalarR,
                                                template_parameters = "Dim<%(ndim)i>::%(VT)s",
                                                pyname = "innerProduct",
                                                cppname = "innerProduct<Dim<%(ndim)i>::%(VT)s>")
""" % {"i" : i,
       "j" : i + 1,
       "VT" : VT,
       "ndim" : ndim})
        i += 2

#-------------------------------------------------------------------------------
# General inner products
#-------------------------------------------------------------------------------
@PYB11template("AType", "BType", "ReturnType")
def innerProduct(A = "const %(AType)s&",
                 B = "const %(BType)s&"):
    "Inner product (%(AType)s.%(BType)s."
    return "%(ReturnType)s"

for AT in ("Vector", "Tensor", "SymTensor", "ThirdRankTensor", "FourthRankTensor", "FifthRankTensor"):
    for BT in ("Vector", "Tensor", "SymTensor", "ThirdRankTensor", "FourthRankTensor", "FifthRankTensor"):
        for ndim in (1, 2, 3):
        exec("""
innerProduct%(i)i = PYB11TemplateFunction(innerProduct,
                                          template_parameters = ("Dim<%(ndim)i>::",
                                          pyname = "innerProduct",
                                          cppname = "innerProduct<Dim<%(ndim)i>")
""" % {"i" : i,
       "AT" : "Dim<
       "VT" : VT,
       "ndim" : ndim})
        i += 1
