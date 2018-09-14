"""
This module provides the fundamental Spheral Geometry types (Vector, Tensor, 
etc.) and associated methods such as products and eigenvalues.
"""

from PYB11Decorators import *
from PYB11STLmethods import *
from PYB11property import *
import types

# Define some useful type collections we're going to be wrapping in this module.
geomtypes = ["Vector", "Tensor", "SymTensor", "ThirdRankTensor", "FourthRankTensor", "FifthRankTensor", "FacetedVolume"]

# We use the preamble for a whole bunch of useful typedefs.
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

            '<vector>']

# STL containers
for gtype in geomtypes + ["Plane"]:
    for suffix in ("1d", "2d", "3d"):
        element = gtype + suffix
        exec('vector_of_%s = PYB11_bind_vector("%s", opaque=True)' % (element, element))
vector_of_Facet2d = PYB11_bind_vector("GeomFacet2d", opaque=True)
vector_of_Facet3d = PYB11_bind_vector("GeomFacet3d", opaque=True)

#-------------------------------------------------------------------------------
# Worker for adding methods to Vector.
#-------------------------------------------------------------------------------
@PYB11ignore
def addVectorMethods(cls, ndim):
    
    me = "Vector%id" % ndim

    # Constructors
    def pyinit0(self):
        "Default constructor"

    def pyinit1(self,
                rhs = "const Vector%id" % ndim):
        "Copy constructor"

    def pyinit2(self,
                x = "double",
                y = ("double", "0.0"),
                z = ("double", "0.0")):
        "Construct with element values."

    # Attributes
    @PYB11static
    @PYB11readonly
    def nDimensions(self):
        "Number of dimensions"

    @PYB11static
    @PYB11readonly
    def numElements(self):
        "Number of elements stored in the type."

    @PYB11static
    @PYB11readonly
    def zero(self):
        "The zero value equivalent."

    @PYB11static
    @PYB11readonly
    def one(self):
        "The unit value equivalent."

    @PYB11cppname("x")
    @PYB11const
    @PYB11ignore
    def getx(self):
        return "double"

    @PYB11cppname("x")
    @PYB11ignore
    def setx(self, val="double"):
        return "void"

    # Properties
    cls.x = PYB11property(getx, setx,
                          doc = "The x coordinate.")

    # Add all the locally defined methods to the cls.
    for _x in [x for x in dir() if type(eval(x)) == types.FunctionType]: 
        exec("cls.%s = %s" % (_x, _x))

#-------------------------------------------------------------------------------
# Vector
#-------------------------------------------------------------------------------
class Vector1d:
    "Spheral Vector (1d)"

addVectorMethods(Vector1d, 1)
