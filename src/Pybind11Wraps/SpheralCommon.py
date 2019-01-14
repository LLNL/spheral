#-------------------------------------------------------------------------------
# Common PYB11 initialization code for all Spheral modules.
#-------------------------------------------------------------------------------
from PYB11Generator import *

PYB11includes = ['"Geometry/Dimension.hh"',
                 '"Geometry/GeomPlane.hh"',
                 "<vector>",
                 "<map>",
                 "<set>",
                 "<string>"]

PYB11preamble = """
typedef Spheral::GeomPlane<Spheral::Dim<1>> Plane1d;
typedef Spheral::Dim<1>::Vector Vector1d;
typedef Spheral::Dim<1>::Tensor Tensor1d;
typedef Spheral::Dim<1>::SymTensor SymTensor1d;
typedef Spheral::Dim<1>::ThirdRankTensor ThirdRankTensor1d;
typedef Spheral::Dim<1>::FourthRankTensor FourthRankTensor1d;
typedef Spheral::Dim<1>::FifthRankTensor FifthRankTensor1d;
typedef Spheral::Dim<1>::FacetedVolume FacetedVolume1d;
typedef Spheral::GeomPlane<Spheral::Dim<2>> Plane2d;
typedef Spheral::Dim<2>::Vector Vector2d;
typedef Spheral::Dim<2>::Tensor Tensor2d;
typedef Spheral::Dim<2>::SymTensor SymTensor2d;
typedef Spheral::Dim<2>::ThirdRankTensor ThirdRankTensor2d;
typedef Spheral::Dim<2>::FourthRankTensor FourthRankTensor2d;
typedef Spheral::Dim<2>::FifthRankTensor FifthRankTensor2d;
typedef Spheral::Dim<2>::FacetedVolume FacetedVolume2d;
typedef Spheral::GeomPlane<Spheral::Dim<3>> Plane3d;
typedef Spheral::Dim<3>::Vector Vector3d;
typedef Spheral::Dim<3>::Tensor Tensor3d;
typedef Spheral::Dim<3>::SymTensor SymTensor3d;
typedef Spheral::Dim<3>::ThirdRankTensor ThirdRankTensor3d;
typedef Spheral::Dim<3>::FourthRankTensor FourthRankTensor3d;
typedef Spheral::Dim<3>::FifthRankTensor FifthRankTensor3d;
typedef Spheral::Dim<3>::FacetedVolume FacetedVolume3d;
"""

# STL containers of common geometric types
from GeometryMOD import geomtypes
for element in geomtypes:
    for ndim in (1, 2, 3):
        exec('''
vector_of_%(mangle)s = PYB11_bind_vector("%(element)s", opaque=True, local=True)
vector_of_vector_of_%(mangle)s = PYB11_bind_vector("std::vector<%(element)s>", opaque=True, local=True)
''' % {"element": "Dim<" + str(ndim) + ">::" + element,
       "mangle" : element + str(ndim) + "d"})
vector_of_Facet2d = PYB11_bind_vector("GeomFacet2d", opaque=True, local=True)
vector_of_Facet3d = PYB11_bind_vector("GeomFacet3d", opaque=True, local=True)
vector_of_Plane1d = PYB11_bind_vector("GeomPlane<Dim<1>>", opaque=True, local=True)
vector_of_Plane2d = PYB11_bind_vector("GeomPlane<Dim<2>>", opaque=True, local=True)
vector_of_Plane3d = PYB11_bind_vector("GeomPlane<Dim<3>>", opaque=True, local=True)
