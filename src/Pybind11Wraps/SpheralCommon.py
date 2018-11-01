#-------------------------------------------------------------------------------
# Common PYB11 initialization code for all Spheral modules.
#-------------------------------------------------------------------------------
from CXXTypesMOD import *

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

PYBIND11_MAKE_OPAQUE(std::vector<Dim<1>::FacetedVolume>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<2>::FacetedVolume>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<3>::FacetedVolume>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<1>::FifthRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<2>::FifthRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<3>::FifthRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<1>::FourthRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<2>::FourthRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<3>::FourthRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<GeomPlane<Dim<1>>>)
PYBIND11_MAKE_OPAQUE(std::vector<GeomPlane<Dim<2>>>)
PYBIND11_MAKE_OPAQUE(std::vector<GeomPlane<Dim<3>>>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<1>::SymTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<2>::SymTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<3>::SymTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<1>::Tensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<2>::Tensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<3>::Tensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<1>::ThirdRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<2>::ThirdRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<3>::ThirdRankTensor>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<1>::Vector>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<2>::Vector>)
PYBIND11_MAKE_OPAQUE(std::vector<Dim<3>::Vector>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<1>::FacetedVolume>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<2>::FacetedVolume>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<3>::FacetedVolume>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<1>::FifthRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<2>::FifthRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<3>::FifthRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<1>::FourthRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<2>::FourthRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<3>::FourthRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<1>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<2>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<3>::SymTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<1>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<2>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<3>::Tensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<1>::ThirdRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<2>::ThirdRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<3>::ThirdRankTensor>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<1>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<2>::Vector>>)
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Dim<3>::Vector>>)
"""
