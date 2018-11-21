#-------------------------------------------------------------------------------
# Common PYB11 initialization code for all Spheral modules.
#-------------------------------------------------------------------------------
from PYB11Generator import *

includes = ['"Geometry/Dimension.hh"',
            "<vector>",
            "<map>",
            "<set>",
            "<string>"]

preamble = """
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

// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FacetedVolume>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FacetedVolume>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FacetedVolume>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FifthRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FifthRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FifthRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FourthRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FourthRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FourthRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<GeomPlane<Spheral::Dim<1>>>)
// PYBIND11_MAKE_OPAQUE(std::vector<GeomPlane<Spheral::Dim<2>>>)
// PYBIND11_MAKE_OPAQUE(std::vector<GeomPlane<Spheral::Dim<3>>>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::SymTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::SymTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::SymTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Tensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::Tensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::Tensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::ThirdRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::ThirdRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::ThirdRankTensor>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Vector>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::Vector>)
// PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::Vector>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::FacetedVolume>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::FacetedVolume>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::FacetedVolume>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::FifthRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::FifthRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::FifthRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::FourthRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::FourthRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::FourthRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::SymTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::SymTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::SymTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::Tensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::Tensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::Tensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::ThirdRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::ThirdRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::ThirdRankTensor>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<1>::Vector>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<2>::Vector>>)
// PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Dim<3>::Vector>>)
"""

# STL containers of primitives

# std::vector
vector_of_char     = PYB11_bind_vector("char", opaque=True, local=True)
vector_of_unsigned = PYB11_bind_vector("unsigned", opaque=True, local=True)
vector_of_ULL      = PYB11_bind_vector("uint64_t", opaque=True, local=True)
vector_of_int      = PYB11_bind_vector("int", opaque=True, local=True)
vector_of_float    = PYB11_bind_vector("float", opaque=True, local=True)
vector_of_double   = PYB11_bind_vector("double", opaque=True, local=True)
vector_of_string   = PYB11_bind_vector("std::string", opaque=True, local=True)

# std::vector<std::vector>
vector_of_vector_of_char     = PYB11_bind_vector("std::vector<char>", opaque=True, local=True)
vector_of_vector_of_unsigned = PYB11_bind_vector("std::vector<unsigned>", opaque=True, local=True)
vector_of_vector_of_ULL      = PYB11_bind_vector("std::vector<uint64_t>", opaque=True, local=True)
vector_of_vector_of_int      = PYB11_bind_vector("std::vector<int>", opaque=True, local=True)
vector_of_vector_of_float    = PYB11_bind_vector("std::vector<float>", opaque=True, local=True)
vector_of_vector_of_double   = PYB11_bind_vector("std::vector<double>", opaque=True, local=True)
vector_of_vector_of_string   = PYB11_bind_vector("std::vector<std::string>", opaque=True, local=True)

# std::vector<pair<>>
vector_of_pair_double_double     = PYB11_bind_vector("std::pair<double, double>", opaque=True, local=True)
vector_of_pair_double_string     = PYB11_bind_vector("std::pair<double, std::string>", opaque=True, local=True)
vector_of_pair_unsigned_unsigned = PYB11_bind_vector("std::pair<unsigned, unsigned>", opaque=True, local=True)
vector_of_pair_ULL_ULL           = PYB11_bind_vector("std::pair<uint64_t, uint64_t>", opaque=True, local=True)
vector_of_pair_string_string     = PYB11_bind_vector("std::pair<std::string, std::string>", opaque=True, local=True)

# std::map
map_string_double = PYB11_bind_map("std::string", "double", opaque=True, local=True)
map_int_string    = PYB11_bind_map("int", "std::string", opaque=True, local=True)

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
