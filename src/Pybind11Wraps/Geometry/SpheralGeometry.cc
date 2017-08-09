// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>

#include "Geometry/Dimension.hh"
#include "Geometry/GeomVector.hh"
#include "Geometry/Geom3Vector.hh"
#include "Geometry/GeomTensor.hh"
#include "Geometry/GeomSymmetricTensor.hh"
#include "Geometry/GeomThirdRankTensor.hh"
#include "Geometry/GeomFourthRankTensor.hh"
#include "Geometry/GeomFifthRankTensor.hh"
#include "Geometry/EigenStruct.hh"
#include "Geometry/computeEigenValues.hh"
#include "Geometry/GeomPlane.hh"
#include "Geometry/GeomPolygon.hh"
#include "Geometry/GeomPolyhedron.hh"
#include "Geometry/GeomFacet2d.hh"
#include "Geometry/GeomFacet3d.hh"
#include "Geometry/invertRankNTensor.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"
#include "Geometry/innerDoubleProduct.hh"
#include "Geometry/aggregateFacetedVolumes.hh"
#include "Field/Field.hh"
#include "Utilities/DataTypeTraits.hh"

namespace py = pybind11;
using namespace pybind11::literals;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Vector>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Tensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::SymTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::ThirdRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FourthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FifthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPlane<Spheral::Dim<1>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Box1d>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Box1d>>);

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::Vector>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::Tensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::SymTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::ThirdRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FourthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<2>::FifthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPlane<Spheral::Dim<2>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomFacet2d>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPolygon>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::GeomPolygon>>);

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::Vector>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::Tensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::SymTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::ThirdRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FourthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<3>::FifthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPlane<Spheral::Dim<3>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomFacet3d>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPolyhedron>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::GeomPolyhedron>>);

namespace {
// Define a local trait struct to find the name.
template<typename TT> struct TensorTraits;

template<> struct TensorTraits<Spheral::GeomTensor<1>> { static std::string name() { return "Tensor"; } };
template<> struct TensorTraits<Spheral::GeomTensor<2>> { static std::string name() { return "Tensor"; } };
template<> struct TensorTraits<Spheral::GeomTensor<3>> { static std::string name() { return "Tensor"; } };

template<> struct TensorTraits<Spheral::GeomSymmetricTensor<1>> { static std::string name() { return "SymTensor"; } };
template<> struct TensorTraits<Spheral::GeomSymmetricTensor<2>> { static std::string name() { return "SymTensor"; } };
template<> struct TensorTraits<Spheral::GeomSymmetricTensor<3>> { static std::string name() { return "SymTensor"; } };

template<> struct TensorTraits<Spheral::GeomThirdRankTensor<1>> { static std::string name() { return "ThirdRankTensor"; } };
template<> struct TensorTraits<Spheral::GeomThirdRankTensor<2>> { static std::string name() { return "ThirdRankTensor"; } };
template<> struct TensorTraits<Spheral::GeomThirdRankTensor<3>> { static std::string name() { return "ThirdRankTensor"; } };

template<> struct TensorTraits<Spheral::GeomFourthRankTensor<1>> { static std::string name() { return "FourthRankTensor"; } };
template<> struct TensorTraits<Spheral::GeomFourthRankTensor<2>> { static std::string name() { return "FourthRankTensor"; } };
template<> struct TensorTraits<Spheral::GeomFourthRankTensor<3>> { static std::string name() { return "FourthRankTensor"; } };

template<> struct TensorTraits<Spheral::GeomFifthRankTensor<1>> { static std::string name() { return "FifthRankTensor"; } };
template<> struct TensorTraits<Spheral::GeomFifthRankTensor<2>> { static std::string name() { return "FifthRankTensor"; } };
template<> struct TensorTraits<Spheral::GeomFifthRankTensor<3>> { static std::string name() { return "FifthRankTensor"; } };

// Also a common method of printing tensors.
template<typename TT>
std::string
tensor_string_repr(const TT& t) {
  std::string result = TensorTraits<TT>::name() + std::to_string(TT::nDimensions) + "d(";
  for (auto val: t) result += (" " + std::to_string(val) + " ");
  result += ")";
  return result;
}
  
//------------------------------------------------------------------------------
// The methods common to Tensor and SymmetricTensor (both rank 2)
//------------------------------------------------------------------------------
template<typename Dimension,
         typename TensorPB11,
         typename OtherTensorPB11,
         typename Tensor,
         typename OtherTensor>
void
tensorBindings(TensorPB11& tensorPB11, OtherTensorPB11& otherPB11,
               Tensor tensor, OtherTensor otherTensor) {

  // Define types.
  typedef typename Dimension::Vector Vector;

  tensorPB11
    // Instance attributes.
    .def_readonly_static("nDimensions", &Tensor::nDimensions)
    .def_readonly_static("numElements", &Tensor::numElements)
    .def_readonly_static("zero", &Tensor::zero)
    .def_readonly_static("one", &Tensor::one)
  
    // Constructors.
    .def(py::init<>())
    .def(py::init<const Tensor&>())
    .def(py::init<const OtherTensor&>())
  
    // Components.
    .def_property("xx", (double (Tensor::*)() const) &Tensor::xx, (void (Tensor::*)(const double)) &Tensor::xx)
    .def_property("xy", (double (Tensor::*)() const) &Tensor::xy, (void (Tensor::*)(const double)) &Tensor::xy)
    .def_property("xz", (double (Tensor::*)() const) &Tensor::xz, (void (Tensor::*)(const double)) &Tensor::xz)
    .def_property("yx", (double (Tensor::*)() const) &Tensor::yx, (void (Tensor::*)(const double)) &Tensor::yx)
    .def_property("yy", (double (Tensor::*)() const) &Tensor::yy, (void (Tensor::*)(const double)) &Tensor::yy)
    .def_property("yz", (double (Tensor::*)() const) &Tensor::yz, (void (Tensor::*)(const double)) &Tensor::yz)
    .def_property("zx", (double (Tensor::*)() const) &Tensor::zx, (void (Tensor::*)(const double)) &Tensor::zx)
    .def_property("zy", (double (Tensor::*)() const) &Tensor::zy, (void (Tensor::*)(const double)) &Tensor::zy)
    .def_property("zz", (double (Tensor::*)() const) &Tensor::zz, (void (Tensor::*)(const double)) &Tensor::zz)
  
    // Add sequence methods.
    .def("__getitem__", [](const Tensor &s, size_t i) {
        if (i >= Tensor::numElements) throw py::index_error();
        return s[i];
      })
    .def("__setitem__", [](Tensor &s, size_t i, float v) {
        if (i >= Tensor::numElements) throw py::index_error();
        s[i] = v;
      })
    .def("__len__", [](const Tensor& self) { return Tensor::numElements; })

    // Optional sequence protocol operations
    .def("__iter__", [](const Tensor &s) { return py::make_iterator(s.begin(), s.end()); },
         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("__contains__", [](const Tensor& self, const double val) { for (const auto elem: self) { if (elem == val) return true; } return false; })

    // We like in index with operator()
    .def("__call__", (double (Tensor::*)(const typename Tensor::size_type,
                                         const typename Tensor::size_type) const) &Tensor::operator(), py::is_operator())

    // Operators.
    .def(-py::self)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self * py::self)
    .def(py::self += py::self)
    .def(py::self -= py::self)

    .def(py::self + OtherTensor())
    .def(py::self - OtherTensor())
    .def(py::self * OtherTensor())

    .def(py::self * float())
    .def(float() * py::self)
    .def(py::self / float())
    .def(py::self *= float())
    .def(py::self /= float())
  
    // Comparisons.
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self)
    .def(py::self >  py::self)
    .def(py::self <= py::self)
    .def(py::self >= py::self)
  
    // Misc methods.
    .def("getRow",&Tensor::getRow)
    .def("getColumn",&Tensor::getColumn)
    .def("Zero", &Tensor::Zero)
    .def("Identity", &Tensor::Identity)
    .def("Symmetric", &Tensor::Symmetric)
    .def("SkewSymmetric", &Tensor::SkewSymmetric)
    .def("Transpose", &Tensor::Transpose)
    .def("Inverse", &Tensor::Inverse)
    .def("diagonalElements", &Tensor::diagonalElements)
    .def("Trace", &Tensor::Trace)
    .def("Determinant", &Tensor::Determinant)
    .def("dot", (Vector (Tensor::*)(const Vector&) const) &Tensor::dot)
    .def("doubledot", (double (Tensor::*)(const Tensor&) const) &Tensor::doubledot)
    .def("doubledot", (double (Tensor::*)(const OtherTensor&) const) &Tensor::doubledot)
    .def("selfDoubledot", &Tensor::selfDoubledot)
    .def("square", &Tensor::square)
    .def("squareElements", &Tensor::squareElements)
    .def("eigenValues", &Tensor::eigenValues)
    .def("rotationalTransform", &Tensor::rotationalTransform)
    .def("maxAbsElement", &Tensor::maxAbsElement)

    // A nicer print.
    .def("__str__", [](const Tensor& self) { return tensor_string_repr(self); })
    .def("__repr__", [](const Tensor& self) { return tensor_string_repr(self); })
    ;

    // Dimension dependent constructors.
    if (Tensor::nDimensions == 1) {
      tensorPB11.def(py::init<double>(), "xx"_a=0.0);
    } else if (Tensor::nDimensions == 2) {
      tensorPB11.def(py::init<double, double, double, double>(),
                     "xx"_a=0.0, "xy"_a=0.0,
                     "yx"_a=0.0, "yy"_a=0.0);
    } else {
      CHECK(Tensor::nDimensions == 3);
      tensorPB11.def(py::init<double, double, double, double, double, double, double, double, double>(),
                     "xx"_a=0.0, "xy"_a=0.0, "xz"_a=0.0,
                     "xy"_a=0.0, "yy"_a=0.0, "yz"_a=0.0,
                     "zx"_a=0.0, "zy"_a=0.0, "zz"_a=0.0);
    }
}

//------------------------------------------------------------------------------
// The methods common to 3rd, 4th, and 5th rank tensors.
//------------------------------------------------------------------------------
template<typename TensorPB11,
         typename Tensor>
void
rankNTensorBindings(TensorPB11& tensorPB11, Tensor tensor) {

  tensorPB11

    // Instance attributes.
    .def_readonly_static("nrank", &Tensor::nrank)
    .def_readonly_static("nDimensions", &Tensor::nDimensions)
    .def_readonly_static("numElements", &Tensor::numElements)
    .def_readonly_static("zero", &Tensor::zero)
    
    // Constructors.
    .def(py::init<>())
    .def(py::init<double>(), "val"_a)
    .def(py::init<const Tensor&>(), "rhs"_a)

    // Add sequence methods.
    .def("__getitem__", [](const Tensor &s, size_t i) {
        if (i >= Tensor::numElements) throw py::index_error();
        return s[i];
      })
    .def("__setitem__", [](Tensor &s, size_t i, float v) {
        if (i >= Tensor::numElements) throw py::index_error();
        s[i] = v;
      })
    .def("__len__", [](const Tensor& self) { return Tensor::numElements; })

    // Optional sequence protocol operations
    .def("__iter__", [](const Tensor &s) { return py::make_iterator(s.begin(), s.end()); },
         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("__contains__", [](const Tensor& self, const double val) { for (const auto elem: self) { if (elem == val) return true; } return false; })

    // A nicer print.
    .def("__str__", [](const Tensor& self) { return tensor_string_repr(self); })
    .def("__repr__", [](const Tensor& self) { return tensor_string_repr(self); })
    ;
}

//------------------------------------------------------------------------------
// Dimesion dependent bindings.
//------------------------------------------------------------------------------
template<int ndim>
void
geometryBindings(py::module& m, const std::string& suffix) {

  // Define types.
  typedef Spheral::Dim<ndim> Dimension;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef Spheral::EigenStruct<ndim> EigenStructType;
  typedef Spheral::GeomPlane<Spheral::Dim<ndim>> PlaneType;

  // Declare the pybind11 types.
  py::class_<Vector> VectorPB11(m, ("Vector" + suffix).c_str(), py::metaclass());
  py::class_<Tensor> TensorPB11(m, ("Tensor" + suffix).c_str(), py::metaclass());
  py::class_<SymTensor> SymTensorPB11(m, ("SymTensor" + suffix).c_str(), py::metaclass());
  py::class_<ThirdRankTensor> ThirdRankTensorPB11(m, ("ThirdRankTensor" + suffix).c_str(), py::metaclass());
  py::class_<FourthRankTensor> FourthRankTensorPB11(m, ("FourthRankTensor" + suffix).c_str(), py::metaclass());
  py::class_<FifthRankTensor> FifthRankTensorPB11(m, ("FifthRankTensor" + suffix).c_str(), py::metaclass());
  py::class_<EigenStructType> EigenStructPB11(m, ("EigenStruct" + suffix).c_str());
  py::class_<PlaneType> PlanePB11(m, ("Plane" + suffix).c_str());

  //............................................................................
  // The STL containers of Geometry objects.
  py::bind_vector<std::vector<Vector>>(m, "vector_of_Vector" + suffix);
  py::bind_vector<std::vector<Tensor>>(m, "vector_of_Tensor" + suffix);
  py::bind_vector<std::vector<SymTensor>>(m, "vector_of_SymTensor" + suffix);
  py::bind_vector<std::vector<ThirdRankTensor>>(m, "vector_of_ThirdRankTensor" + suffix);
  py::bind_vector<std::vector<FourthRankTensor>>(m, "vector_of_FourthRankTensor" + suffix);
  py::bind_vector<std::vector<FifthRankTensor>>(m, "vector_of_FifthRankTensor" + suffix);
  py::bind_vector<std::vector<PlaneType>>(m, "vector_of_Plane" + suffix);
  py::bind_vector<std::vector<typename Dimension::FacetedVolume>>(m, "vector_of_FacetedVolume" + suffix);

  //............................................................................
  // Vector
  VectorPB11

    // Instance attributes.
    .def_readonly_static("nDimensions", &Vector::nDimensions)
    .def_readonly_static("numElements", &Vector::numElements)
    .def_readonly_static("zero", &Vector::zero)
    .def_readonly_static("one", &Vector::one)
    
    // Constructors.
    .def(py::init<>())
    .def(py::init<const Vector&>())
    .def(py::init<double, double, double>(), "x"_a, "y"_a=0.0, "z"_a=0.0)
    // .def(py::init<double, double, double>(), "x"_a, "y"_a=0.0, "z"_a=0.0)
    
    // x, y, z
    .def_property("x", (double (Vector::*)() const) &Vector::x, (void (Vector::*)(const double)) &Vector::x)
    .def_property("y", (double (Vector::*)() const) &Vector::y, (void (Vector::*)(const double)) &Vector::y)
    .def_property("z", (double (Vector::*)() const) &Vector::z, (void (Vector::*)(const double)) &Vector::z)
    
    // Add sequence methods.
    .def("__getitem__", [](const Vector &s, size_t i) {
        if (i >= Vector::numElements) throw py::index_error();
        return s[i];
      })
    .def("__setitem__", [](Vector &s, size_t i, float v) {
        if (i >= Vector::numElements) throw py::index_error();
        s[i] = v;
      })
    .def("__len__", [](const Vector& self) { return Dimension::nDim; })

    // Optional sequence protocol operations
    .def("__iter__", [](const Vector &s) { return py::make_iterator(s.begin(), s.end()); },
         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)

    // We like in index with operator()
    .def("__call__", (double (Vector::*)(typename Vector::size_type) const) &Vector::operator(), py::is_operator())

    .def("Zero", &Vector::Zero)
    
    // Operators.
    .def(-py::self)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self * py::self)
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self * float())
    .def(float() * py::self)
    .def(py::self / float())
    .def(py::self *= float())
    .def(py::self /= float())
    
    // Comparisons.
    .def("compare", (int (Vector::*)(const Vector&) const)  &Vector::compare)
    .def("compare", (int (Vector::*)(const double) const)  &Vector::compare)
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self)
    .def(py::self >  py::self)
    .def(py::self <= py::self)
    .def(py::self >= py::self)
    
    // Misc methods.
    .def("dot", &Vector::dot)
    .def("cross", &Vector::cross)
    .def("dyad", &Vector::dyad)
    .def("selfdyad", &Vector::selfdyad)
    .def("unitVector", &Vector::unitVector)
    .def("magnitude", &Vector::magnitude)
    .def("magnitude2", &Vector::magnitude2)
    .def("minElement", &Vector::minElement)
    .def("maxElement", &Vector::maxElement)
    .def("maxAbsElement", &Vector::maxAbsElement)
    .def("sumElements", &Vector::sumElements)

    // A nicer print.
    .def("__str__", [](const Vector& self) {
        std::string result = "Vector" + std::to_string(Dimension::nDim) + "d(";
        for (auto val: self) result += (" " + std::to_string(val) + " ");
        result += ")";
        return result;
      }
      )
    .def("__repr__", [](const Vector& self) {
        std::string result = "Vector" + std::to_string(Dimension::nDim) + "d(";
        for (auto val: self) result += (" " + std::to_string(val) + " ");
        result += ")";
        return result;
      }
      )
    ;

  //............................................................................
  // Tensor
  tensorBindings<Dimension>(TensorPB11, SymTensorPB11, Tensor(), SymTensor());
  TensorPB11
    .def("dot", (Tensor (Tensor::*)(const Tensor&) const) &Tensor::dot)
    .def("dot", (Tensor (Tensor::*)(const SymTensor&) const) &Tensor::dot)
    .def(py::self += SymTensor())
    .def(py::self -= SymTensor())
    ;

  //............................................................................
  // SymTensor
  tensorBindings<Dimension>(SymTensorPB11, TensorPB11, SymTensor(), Tensor());
  SymTensorPB11
    .def("cube", &SymTensor::cube)
    .def("sqrt", &SymTensor::sqrt)
    .def("cuberoot", &SymTensor::cuberoot)
    .def("pow", &SymTensor::pow)
    .def("eigenVectors", &SymTensor::eigenVectors)
    ;

  //............................................................................
  // ThirdRankTensor
  rankNTensorBindings(ThirdRankTensorPB11, ThirdRankTensor());
  ThirdRankTensorPB11
    .def("__call__", (double (ThirdRankTensor::*)(const typename ThirdRankTensor::size_type,
                                                  const typename ThirdRankTensor::size_type,
                                                  const typename ThirdRankTensor::size_type) const) &ThirdRankTensor::operator(), py::is_operator())
    .def("__call__", [](ThirdRankTensor& self, const int i, const int j, const int k, const double val){
        VERIFY2(i < ThirdRankTensor::nDimensions and
                j < ThirdRankTensor::nDimensions and
                k < ThirdRankTensor::nDimensions,
                "ThirdRankTensor assignment error : out of bounds index");
        self(i,j,k) = val;
      }
      )
    ;

  //............................................................................
  // FourthRankTensor
  rankNTensorBindings(FourthRankTensorPB11, FourthRankTensor());
  FourthRankTensorPB11
    .def("__call__", (double (FourthRankTensor::*)(const typename FourthRankTensor::size_type,
                                                   const typename FourthRankTensor::size_type,
                                                   const typename FourthRankTensor::size_type,
                                                   const typename FourthRankTensor::size_type) const) &FourthRankTensor::operator(), py::is_operator())
    .def("__call__", [](FourthRankTensor& self, const int i, const int j, const int k, const int m, const double val){
        VERIFY2(i < FourthRankTensor::nDimensions and
                j < FourthRankTensor::nDimensions and
                k < FourthRankTensor::nDimensions and
                m < FourthRankTensor::nDimensions,
                "FourthRankTensor assignment error : out of bounds index");
        self(i,j,k,m) = val;
      }
      )
    ;

  //............................................................................
  // FifthRankTensor
  rankNTensorBindings(FifthRankTensorPB11, FifthRankTensor());
  FifthRankTensorPB11
    .def("__call__", (double (FifthRankTensor::*)(const typename FifthRankTensor::size_type,
                                                  const typename FifthRankTensor::size_type,
                                                  const typename FifthRankTensor::size_type,
                                                  const typename FifthRankTensor::size_type,
                                                  const typename FifthRankTensor::size_type) const) &FifthRankTensor::operator(), py::is_operator())
    .def("__call__", [](FifthRankTensor& self, const int i, const int j, const int k, const int m, const int n, const double val){
        VERIFY2(i < FifthRankTensor::nDimensions and
                j < FifthRankTensor::nDimensions and
                k < FifthRankTensor::nDimensions and
                m < FifthRankTensor::nDimensions and
                n < FifthRankTensor::nDimensions,
                "FifthRankTensor assignment error : out of bounds index");
        self(i,j,k,m,n) = val;
      }
      )
    ;

  //............................................................................
  // EigenStruct
  EigenStructPB11

    .def(py::init<>())
    .def(py::init<const EigenStructType&>(), "rhs"_a)
    .def_readwrite("eigenValues", &EigenStructType::eigenValues)
    .def_readwrite("eigenVectors", &EigenStructType::eigenVectors)
    // A nicer print.
    .def("__str__", [](const EigenStructType& self) {
        std::stringstream ss;
        ss << "EigenStruct" << Dimension::nDim << "d[" << self.eigenValues << ", " << self.eigenVectors << "]";
        return ss.str();
      }
      )
    .def("__repr__", [](const EigenStructType& self) {
        std::stringstream ss;
        ss << "EigenStruct" << Dimension::nDim << "d[" << self.eigenValues << ", " << self.eigenVectors << "]";
        return ss.str();
      }
      )
    ;

  //............................................................................
  // Plane
  PlanePB11
    .def(py::init<>())
    .def(py::init<const PlaneType&>(), "rhs"_a)
    .def(py::init<const Vector&, const Vector&>(), "point"_a, "normal"_a)
    .def(py::init<const std::vector<Vector>&>(), "points"_a)
    .def_property("point",
                  (const Vector& (PlaneType::*)() const) &PlaneType::point, 
                  (void (PlaneType::*)(const Vector&)) &PlaneType::point)
    .def_property("normal",
                  (const Vector& (PlaneType::*)() const) &PlaneType::normal, 
                  (void (PlaneType::*)(const Vector&)) &PlaneType::normal)
    .def("signedDistance", &PlaneType::signedDistance)
    .def("minimumDistance", &PlaneType::signedDistance)
    .def("closestPointOnPlane", &PlaneType::signedDistance)
    .def("parallel", &PlaneType::parallel)
    .def("valid", &PlaneType::valid)
    .def("compare", &PlaneType::compare)
    .def(-py::self)
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self)
    ;

  //............................................................................
  // Geometry functions
  m.def("invertRankNTensor", Spheral::invertRankNTensor<Tensor>);
  m.def("invertRankNTensor", Spheral::invertRankNTensor<SymTensor>);
  m.def("invertRankNTensor", Spheral::invertRankNTensor<FourthRankTensor>);
  m.def("computeEigenValues", Spheral::computeEigenValues<Dimension>);

  // Inner product (double)
  m.def("innerProduct", (Vector (*)(const double&, const Vector&)) Geometry::innerProduct<Vector>);
  m.def("innerProduct", (Vector (*)(const Vector&, const double&)) Geometry::innerProduct<Vector>);
  m.def("innerProduct", (Tensor (*)(const double&, const Tensor&)) Geometry::innerProduct<Tensor>);
  m.def("innerProduct", (Tensor (*)(const Tensor&, const double&)) Geometry::innerProduct<Tensor>);
  m.def("innerProduct", (SymTensor (*)(const double&, const SymTensor&)) Geometry::innerProduct<SymTensor>);
  m.def("innerProduct", (SymTensor (*)(const SymTensor&, const double&)) Geometry::innerProduct<SymTensor>);
  m.def("innerProduct", (ThirdRankTensor (*)(const double&, const ThirdRankTensor&)) Geometry::innerProduct<ThirdRankTensor>);
  m.def("innerProduct", (ThirdRankTensor (*)(const ThirdRankTensor&, const double&)) Geometry::innerProduct<ThirdRankTensor>);
  m.def("innerProduct", (FourthRankTensor (*)(const double&, const FourthRankTensor&)) Geometry::innerProduct<FourthRankTensor>);
  m.def("innerProduct", (FourthRankTensor (*)(const FourthRankTensor&, const double&)) Geometry::innerProduct<FourthRankTensor>);
  m.def("innerProduct", (FifthRankTensor (*)(const double&, const FifthRankTensor&)) Geometry::innerProduct<FifthRankTensor>);
  m.def("innerProduct", (FifthRankTensor (*)(const FifthRankTensor&, const double&)) Geometry::innerProduct<FifthRankTensor>);

  // Inner product (Vector)
  m.def("innerProduct", (double (*)(const Vector&, const Vector&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Vector (*)(const Vector&, const Tensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Vector (*)(const Tensor&, const Vector&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Vector (*)(const Vector&, const SymTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Vector (*)(const SymTensor&, const Vector&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Tensor (*)(const Vector&, const ThirdRankTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Tensor (*)(const ThirdRankTensor&, const Vector&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (ThirdRankTensor (*)(const Vector&, const FourthRankTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (ThirdRankTensor (*)(const FourthRankTensor&, const Vector&)) Geometry::innerProduct<Dimension>);

  // Inner product (Tensor)
  m.def("innerProduct", (Tensor (*)(const Tensor&, const Tensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Tensor (*)(const Tensor&, const SymTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (Tensor (*)(const SymTensor&, const Tensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (ThirdRankTensor (*)(const Tensor&, const ThirdRankTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (ThirdRankTensor (*)(const ThirdRankTensor&, const Tensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (FourthRankTensor (*)(const Tensor&, const FourthRankTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (FourthRankTensor (*)(const FourthRankTensor&, const Tensor&)) Geometry::innerProduct<Dimension>);

  // Inner product (ThirdRankTensor)
  m.def("innerProduct", (FourthRankTensor (*)(const ThirdRankTensor&, const ThirdRankTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (FifthRankTensor (*)(const ThirdRankTensor&, const FourthRankTensor&)) Geometry::innerProduct<Dimension>);
  m.def("innerProduct", (FifthRankTensor (*)(const FourthRankTensor&, const ThirdRankTensor&)) Geometry::innerProduct<Dimension>);

  // Outer product (double)
  m.def("outerProduct", (Vector (*)(const double&, const Vector&)) Geometry::outerProduct<Vector>);
  m.def("outerProduct", (Vector (*)(const Vector&, const double&)) Geometry::outerProduct<Vector>);
  m.def("outerProduct", (Tensor (*)(const double&, const Tensor&)) Geometry::outerProduct<Tensor>);
  m.def("outerProduct", (Tensor (*)(const Tensor&, const double&)) Geometry::outerProduct<Tensor>);
  m.def("outerProduct", (SymTensor (*)(const double&, const SymTensor&)) Geometry::outerProduct<SymTensor>);
  m.def("outerProduct", (SymTensor (*)(const SymTensor&, const double&)) Geometry::outerProduct<SymTensor>);
  m.def("outerProduct", (ThirdRankTensor (*)(const double&, const ThirdRankTensor&)) Geometry::outerProduct<ThirdRankTensor>);
  m.def("outerProduct", (ThirdRankTensor (*)(const ThirdRankTensor&, const double&)) Geometry::outerProduct<ThirdRankTensor>);
  m.def("outerProduct", (FourthRankTensor (*)(const double&, const FourthRankTensor&)) Geometry::outerProduct<FourthRankTensor>);
  m.def("outerProduct", (FourthRankTensor (*)(const FourthRankTensor&, const double&)) Geometry::outerProduct<FourthRankTensor>);
  m.def("outerProduct", (FifthRankTensor (*)(const double&, const FifthRankTensor&)) Geometry::outerProduct<FifthRankTensor>);
  m.def("outerProduct", (FifthRankTensor (*)(const FifthRankTensor&, const double&)) Geometry::outerProduct<FifthRankTensor>);

  // Outer product (Vector)
  m.def("outerProduct", (Tensor (*)(const Vector&, const Vector&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (ThirdRankTensor (*)(const Vector&, const Tensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (ThirdRankTensor (*)(const Tensor&, const Vector&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (ThirdRankTensor (*)(const Vector&, const SymTensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (ThirdRankTensor (*)(const SymTensor&, const Vector&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FourthRankTensor (*)(const Vector&, const ThirdRankTensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FourthRankTensor (*)(const ThirdRankTensor&, const Vector&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FifthRankTensor (*)(const Vector&, const FourthRankTensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FifthRankTensor (*)(const FourthRankTensor&, const Vector&)) Geometry::outerProduct<Dimension>);

  // Outer product (Tensor)
  m.def("outerProduct", (FourthRankTensor (*)(const Tensor&, const Tensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FourthRankTensor (*)(const Tensor&, const SymTensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FourthRankTensor (*)(const SymTensor&, const Tensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FifthRankTensor (*)(const Tensor&, const ThirdRankTensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FifthRankTensor (*)(const ThirdRankTensor&, const Tensor&)) Geometry::outerProduct<Dimension>);

  // Outer product (SymTensor)
  m.def("outerProduct", (FourthRankTensor (*)(const SymTensor&, const SymTensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FifthRankTensor (*)(const SymTensor&, const ThirdRankTensor&)) Geometry::outerProduct<Dimension>);
  m.def("outerProduct", (FifthRankTensor (*)(const ThirdRankTensor&, const SymTensor&)) Geometry::outerProduct<Dimension>);

  // Inner double product (Tensor)
  m.def("innerDoubleProduct", (double (*)(const Tensor&, const Tensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (double (*)(const Tensor&, const SymTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (double (*)(const SymTensor&, const Tensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Vector (*)(const Tensor&, const ThirdRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Vector (*)(const ThirdRankTensor&, const Tensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Tensor (*)(const Tensor&, const FourthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Tensor (*)(const FourthRankTensor&, const Tensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (ThirdRankTensor (*)(const Tensor&, const FifthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (ThirdRankTensor (*)(const FifthRankTensor&, const Tensor&)) Geometry::innerDoubleProduct<Dimension>);

  // Inner double product (SymTensor)
  m.def("innerDoubleProduct", (double (*)(const SymTensor&, const SymTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Vector (*)(const SymTensor&, const ThirdRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Vector (*)(const ThirdRankTensor&, const SymTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Tensor (*)(const SymTensor&, const FourthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (Tensor (*)(const FourthRankTensor&, const SymTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (ThirdRankTensor (*)(const SymTensor&, const FifthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (ThirdRankTensor (*)(const FifthRankTensor&, const SymTensor&)) Geometry::innerDoubleProduct<Dimension>);

  // Inner double product (ThirdRankTensor)
  m.def("innerDoubleProduct", (Tensor (*)(const ThirdRankTensor&, const ThirdRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (ThirdRankTensor (*)(const ThirdRankTensor&, const FourthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (ThirdRankTensor (*)(const FourthRankTensor&, const ThirdRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (FourthRankTensor (*)(const ThirdRankTensor&, const FifthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (FourthRankTensor (*)(const FifthRankTensor&, const ThirdRankTensor&)) Geometry::innerDoubleProduct<Dimension>);

  // Inner double product (FourthRankTensor)
  m.def("innerDoubleProduct", (FourthRankTensor (*)(const FourthRankTensor&, const FourthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (FifthRankTensor (*)(const FourthRankTensor&, const FifthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);
  m.def("innerDoubleProduct", (FifthRankTensor (*)(const FifthRankTensor&, const FourthRankTensor&)) Geometry::innerDoubleProduct<Dimension>);

}
}

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralGeometry) {
  py::module m("SpheralGeometry", "Spheral Geometry module types.");

  //............................................................................
  // Generic dimension objects.
  geometryBindings<1>(m, "1d");
  geometryBindings<2>(m, "2d");
  geometryBindings<3>(m, "3d");

  //............................................................................
  // Box1d
  typedef Spheral::Box1d Box;
  typedef Spheral::GeomVector<1> Vector1d;
  py::class_<Box>(m, "Box1d")

    // Constructors
    .def(py::init<>())
    .def(py::init<const std::vector<Vector1d>&>(), "points"_a)
    .def(py::init<const std::vector<Vector1d>&, const std::vector<std::vector<unsigned>>>(), "points"_a, "facetIndices"_a)
    .def(py::init<const Vector1d&, const double>(), "center"_a, "extent"_a)
    .def(py::init<const Box&>(), "rhs"_a)

    // Attributes
    .def_property("center", (const Vector1d& (Box::*)() const) &Box::center, (void (Box::*)(const Vector1d&)) &Box::center)
    .def_property("extent", (double (Box::*)() const) &Box::extent, (void (Box::*)(const double)) &Box::extent)
    .def_property_readonly("xmin", &Box::xmin)
    .def_property_readonly("xmax", &Box::xmax)
    .def_property_readonly("volume", &Box::volume)

    // Methods
    .def("contains", &Box::contains)
    .def("convexContains", &Box::convexContains)
    .def("intersect", (bool (Box::*)(const Box&) const) &Box::intersect)
    .def("intersect", (bool (Box::*)(const std::pair<Vector1d, Vector1d>&) const) &Box::intersect)
    .def("convexIntersect", &Box::convexIntersect)
    .def("distance", &Box::distance)
    .def("closestPoint", &Box::closestPoint)
    .def("vertices", &Box::vertices)

    // Operators
    .def(py::self += Vector1d())
    .def(py::self -= Vector1d())
    .def(py::self + Vector1d())
    .def(py::self - Vector1d())

    .def(py::self *= float())
    .def(py::self /= float())
    .def(py::self * float())
    .def(py::self / float())
    
    .def(py::self == py::self)
    .def(py::self != py::self)
    ;

  //............................................................................
  // Polygon
  typedef Spheral::GeomPolygon Polygon;
  typedef Spheral::GeomVector<2> Vector2d;
  py::class_<Polygon>(m, "Polygon")

    // Constructors
    .def(py::init<>())
    .def(py::init<const std::vector<Vector2d>&>(), "points"_a)
    .def(py::init<const std::vector<Vector2d>&, const std::vector<std::vector<unsigned>>>(), "points"_a, "facetIndices"_a)
    .def(py::init<const Polygon&>(), "rhs"_a)

    // Attributes
    .def_property_readonly("xmin", &Polygon::xmin)
    .def_property_readonly("xmax", &Polygon::xmax)
    .def_property_readonly("edges", &Polygon::edges)
    .def_property_readonly("facetVertices", &Polygon::facetVertices)
    .def_property_readonly("volume", &Polygon::volume)

    // Methods
    .def("contains", &Polygon::contains)
    .def("convexContains", &Polygon::convexContains)
    .def("intersect", (bool (Polygon::*)(const Polygon&) const) &Polygon::intersect)
    .def("intersect", (bool (Polygon::*)(const std::pair<Vector2d, Vector2d>&) const) &Polygon::intersect)
    .def("intersect", (std::vector<Vector2d> (Polygon::*)(const Vector2d&, const Vector2d&) const) &Polygon::intersect)
    .def("convexIntersect", &Polygon::convexIntersect)
    .def("centroid", &Polygon::centroid)
    .def("vertices", &Polygon::vertices)
    .def("facets", &Polygon::facets)
    .def("vertexUnitNorms", &Polygon::vertexUnitNorms)
    .def("vertexFacetConnectivity", &Polygon::vertexFacetConnectivity)
    .def("facetFacetConnectivity", &Polygon::facetFacetConnectivity)
    .def("closestFacet", &Polygon::closestFacet)
    .def("distance", &Polygon::distance)
    .def("closestPoint", &Polygon::closestPoint)
    .def("reconstruct", &Polygon::reconstruct)
    .def("convex", &Polygon::convex)
    .def("setBoundingBox", &Polygon::setBoundingBox)

    // Operators
    .def(py::self += Vector2d())
    .def(py::self -= Vector2d())
    .def(py::self + Vector2d())
    .def(py::self - Vector2d())

    .def(py::self *= float())
    .def(py::self /= float())
    .def(py::self * float())
    .def(py::self / float())
    
    .def(py::self == py::self)
    .def(py::self != py::self)
    ;

  //............................................................................
  // Polyhedron
  typedef Spheral::GeomPolyhedron Polyhedron;
  typedef Spheral::GeomVector<3> Vector3d;
  py::class_<Polyhedron>(m, "Polyhedron")

    // Constructors
    .def(py::init<>())
    .def(py::init<const std::vector<Vector3d>&>(), "points"_a)
    .def(py::init<const std::vector<Vector3d>&, const std::vector<std::vector<unsigned>>>(), "points"_a, "facetIndices"_a)
    .def(py::init<const Polyhedron&>(), "rhs"_a)

    // Attributes
    .def_property_readonly("xmin", &Polyhedron::xmin)
    .def_property_readonly("xmax", &Polyhedron::xmax)
    .def_property_readonly("edges", &Polyhedron::edges)
    .def_property_readonly("facetVertices", &Polyhedron::facetVertices)
    .def_property_readonly("facetNormals", &Polyhedron::facetNormals)
    .def_property_readonly("volume", &Polyhedron::volume)

    // Methods
    .def("contains", &Polyhedron::contains)
    .def("convexContains", &Polyhedron::convexContains)
    .def("intersect", (bool (Polyhedron::*)(const Polyhedron&) const) &Polyhedron::intersect)
    .def("intersect", (bool (Polyhedron::*)(const std::pair<Vector3d, Vector3d>&) const) &Polyhedron::intersect)
    .def("convexIntersect", &Polyhedron::convexIntersect)
    .def("centroid", &Polyhedron::centroid)
    .def("vertices", &Polyhedron::vertices)
    .def("facets", &Polyhedron::facets)
    .def("vertexUnitNorms", &Polyhedron::vertexUnitNorms)
    .def("vertexFacetConnectivity", &Polyhedron::vertexFacetConnectivity)
    .def("facetFacetConnectivity", &Polyhedron::facetFacetConnectivity)
    .def("closestFacet", &Polyhedron::closestFacet)
    .def("distance", &Polyhedron::distance)
    .def("closestPoint", &Polyhedron::closestPoint)
    .def("reconstruct", &Polyhedron::reconstruct)
    .def("convex", &Polyhedron::convex)
    .def("setBoundingBox", &Polyhedron::setBoundingBox)

    // Operators
    .def(py::self += Vector3d())
    .def(py::self -= Vector3d())
    .def(py::self + Vector3d())
    .def(py::self - Vector3d())

    .def(py::self *= float())
    .def(py::self /= float())
    .def(py::self * float())
    .def(py::self / float())
    
    .def(py::self == py::self)
    .def(py::self != py::self)
    ;

  //............................................................................
  // aggregateFacetedVolumes
  m.def("aggregateFacetedVolumes", &Spheral::aggregateFacetedVolumes<Spheral::Dim<2>>);
  m.def("aggregateFacetedVolumes", &Spheral::aggregateFacetedVolumes<Spheral::Dim<3>>);

  //............................................................................
  // The STL containers of 2D and 3D facets.
  py::bind_vector<std::vector<Spheral::GeomFacet2d>>(m, "vector_of_Facet2d");
  py::bind_vector<std::vector<Spheral::GeomFacet3d>>(m, "vector_of_Facet3d");

  return m.ptr();
}
