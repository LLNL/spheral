"""
This module provides the fundamental Spheral Geometry types (Vector, Tensor, 
etc.) and associated methods such as products and eigenvalues.
"""

from PYB11Decorators import *
from PYB11STLmethods import *
from PYB11property import *
from PYB11class import *
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

            '<vector>']

# STL containers
for gtype in geomtypes + ["Plane"]:
    for suffix in ("1d", "2d", "3d"):
        element = gtype + suffix
        exec('vector_of_%s = PYB11_bind_vector("%s", opaque=True)' % (element, element))
vector_of_Facet2d = PYB11_bind_vector("GeomFacet2d", opaque=True)
vector_of_Facet3d = PYB11_bind_vector("GeomFacet3d", opaque=True)

# Get the objects wrapped in split off files.
from Vector import Vector1d, Vector2d, Vector3d
from Tensor import Tensor1d, Tensor2d, Tensor3d
from SymTensor import SymTensor1d, SymTensor2d, SymTensor3d

