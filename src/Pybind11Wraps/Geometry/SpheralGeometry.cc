#include <vector>
#include <map>
#include <set>
#include <string>
#include <sstream>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

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
#include "Utilities/DataTypeTraits.hh"

namespace py = pybind11;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
#ifdef SPHERAL1D
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Vector>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::Tensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::SymTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::ThirdRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FourthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Dim<1>::FifthRankTensor>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::GeomPlane<Spheral::Dim<1>>>);
PYBIND11_MAKE_OPAQUE(std::vector<Spheral::Box1d>);
PYBIND11_MAKE_OPAQUE(std::vector<std::vector<Spheral::Box1d>>);
#endif

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
#ifdef SPHERAL2D
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
#endif

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
#ifdef SPHERAL3D
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
#endif

//------------------------------------------------------------------------------
// Dimesion dependent bindings.
//------------------------------------------------------------------------------
namespace {
template<typename Dimension>
void
geometryBindings(const py::module& m, const std::string& suffix) {

  char namebuffer[128];

  // Define types.
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;

  py::class_<Vector> VectorPB11(m, ("Vector" + suffix).c_str(), py::metaclass());
  py::class_<Tensor> TensorPB11(m, ("Tensor" + suffix).c_str(), py::metaclass());
  py::class_<SymTensor> SymTensorPB11(m, ("SymTensor" + suffix).c_str(), py::metaclass());
  py::class_<ThirdRankTensor> ThirdRankTensorPB11(m, ("ThirdRankTensor" + suffix).c_str(), py::metaclass());
  py::class_<FourthRankTensor> FourthRankTensorPB11(m, ("FourthRankTensor" + suffix).c_str(), py::metaclass());
  py::class_<FifthRankTensor> FifthRankTensorPB11(m, ("FifthRankTensor" + suffix).c_str(), py::metaclass());

  // Vector
  VectorPB11

    // Instance attributes.
    .def_readonly_static("nDimensions", &Vector::nDimensions)
    .def_readonly_static("numElements", &Vector::numElements)
    .def_readonly_static("zero", &Vector::zero)
    .def_readonly_static("one", &Vector::one)
    
    // Constructors.
    .def(py::init<>())
    .def(py::init<double, double, double>(), py::arg("x"), py::arg("y")=0.0, py::arg("z")=0.0)
    // .def(py::init<double, double, double>(), "x"_a, "y"_a=0.0, "z"_a=0.0)
    
    // x, y, z
    .def_property("x", (double (Vector::*)() const) &Vector::x, (void (Vector::*)(const double)) &Vector::x)
    .def_property("y", (double (Vector::*)() const) &Vector::y, (void (Vector::*)(const double)) &Vector::y)
    .def_property("z", (double (Vector::*)() const) &Vector::z, (void (Vector::*)(const double)) &Vector::z)
    
    // Add sequence methods.
    .def("__getitem__", [](const Vector &s, size_t i) {
        if (i >= Dimension::nDim) throw py::index_error();
        return s[i];
      })
    .def("__setitem__", [](Vector &s, size_t i, float v) {
        if (i >= Dimension::nDim) throw py::index_error();
        s[i] = v;
      })
    .def("__len__", []() { return Dimension::nDim; })

    // Optional sequence protocol operations
    .def("__iter__", [](const Vector &s) { return py::make_iterator(s.begin(), s.end()); },
         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)

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
}
}

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralGeometry) {
  py::module m("SpheralGeometry", "Spheral Geometry module types.");

#ifdef SPHERAL1D
  geometryBindings<Spheral::Dim<1>>(m, "1d");
#endif

  return m.ptr();
}
