#include <vector>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "Geometry/Dimension.hh"
#include "Neighbor/GridCellIndex.hh"
#include "Neighbor/GridCellPlane.hh"
#include "Neighbor/Neighbor.hh"
#include "Neighbor/NestedGridNeighbor.hh"
#include "Neighbor/TreeNeighbor.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::NeighborSpace;

namespace {  // anonymous

//------------------------------------------------------------------------------
// Common methods to Neighbor objects.
//------------------------------------------------------------------------------
// template<typename Dimension, typename Obj, typename PB11Obj>
// void genericNeighborBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

//   typedef typename Dimension::Vector Vector;
//   typedef typename Dimension::SymTensor SymTensor;

//   obj

//     // Methods
//     .def("__call__", (double (Obj::*)(double, const double&) const) &Obj::operator(), py::is_operator())
//     .def("__call__", (double (Obj::*)(const Vector&, const double&) const) &Obj::operator(), py::is_operator())
//     .def("__call__", (double (Obj::*)(double, const SymTensor&) const) &Obj::operator(), py::is_operator())
//     .def("__call__", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::operator(), py::is_operator())

//     .def("grad", (double (Obj::*)(double, const double&) const) &Obj::grad)
//     .def("grad", (double (Obj::*)(const Vector&, const double&) const) &Obj::grad)
//     .def("grad", (double (Obj::*)(double, const SymTensor&) const) &Obj::grad)
//     .def("grad", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::grad)

//     .def("grad2", (double (Obj::*)(double, const double&) const) &Obj::grad2)
//     .def("grad2", (double (Obj::*)(const Vector&, const double&) const) &Obj::grad2)
//     .def("grad2", (double (Obj::*)(double, const SymTensor&) const) &Obj::grad2)
//     .def("grad2", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::grad2)

//     .def("gradh", (double (Obj::*)(double, const double&) const) &Obj::gradh)
//     .def("gradh", (double (Obj::*)(const Vector&, const double&) const) &Obj::gradh)
//     .def("gradh", (double (Obj::*)(double, const SymTensor&) const) &Obj::gradh)
//     .def("gradh", (double (Obj::*)(const Vector&, const SymTensor&) const) &Obj::gradh)

//     .def("kernelValue", &Obj::kernelValue)
//     .def("gradValue", &Obj::gradValue)
//     .def("grad2Value", &Obj::grad2Value)
//     .def("gradhValue", &Obj::gradhValue)
//     .def("valid", &Obj::valid)

//     // These methods are in protected scope
//     // .def("setVolumeNormalization", &Obj::setVolumeNormalization)
//     // .def("setKernelExtent", &Obj::setKernelExtent)
//     // .def("setInflectionPoint", &Obj::setInflectionPoint)

//     // Attributes
//     .def_property_readonly("volumeNormalization", &Obj::volumeNormalization)
//     .def_property_readonly("kernelExtent", &Obj::kernelExtent)
//     .def_property_readonly("inflectionPoint", &Obj::inflectionPoint)
//     ;
// }

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  //............................................................................
  // GridCellIndex
  typedef GridCellIndex<Dimension> GCI;
  py::class_<GCI>(m, ("GridCellIndex" + suffix).c_str())

    // Constructors
    .def(py::init<>())
    .def(py::init<int>(), "xIndex"_a)
    .def(py::init<int, int>(), "xIndex"_a, "yIndex"_a)
    .def(py::init<int, int, int>(), "xIndex"_a, "yIndex"_a, "zIndex"_a)
    .def(py::init<const GCI&>(), "rhs"_a)

    // Methods
    .def("setIndices", (void (GCI::*)(int)) &GCI::setIndices, "xIndex"_a)
    .def("setIndices", (void (GCI::*)(int, int)) &GCI::setIndices, "xIndex"_a, "yIndex"_a)
    .def("setIndices", (void (GCI::*)(int, int, int)) &GCI::setIndices, "xIndex"_a, "yIndex"_a, "zIndex"_a)
    .def("dot", &GCI::dot)
    .def("compare", &GCI::compare)
    .def("inRange", &GCI::inRange)
    .def("magnitude", &GCI::magnitude)
    .def("minElement", &GCI::minElement)
    .def("maxElement", &GCI::maxElement)
    .def("sumElements", &GCI::sumElements)
    .def("productElements", &GCI::productElements)
    .def("indexMin", &GCI::indexMin)
    .def("indexMax", &GCI::indexMax)

    // Operators
    .def("__call__", (int (GCI::*)(int) const) &GCI::operator(), "index"_a, py::is_operator())
    .def(-py::self)
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self + int())
    .def(py::self - int())
    .def(py::self * int())
    .def(py::self / int())

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self <  py::self)
    .def(py::self >  py::self)
    .def(py::self <= py::self)
    .def(py::self >= py::self)

    // Attributes
    .def_property("xIndex", (int (GCI::*)() const) &GCI::xIndex, (void (GCI::*)(int)) &GCI::xIndex)
    .def_property("yIndex", (int (GCI::*)() const) &GCI::yIndex, (void (GCI::*)(int)) &GCI::yIndex)
    .def_property("zIndex", (int (GCI::*)() const) &GCI::zIndex, (void (GCI::*)(int)) &GCI::zIndex)
    ;

  //............................................................................
  // GridCellPlane
  typedef GridCellPlane<Dimension> GCP;
  py::class_<GCP>(m, ("GridCellPlane" + suffix).c_str())

    // Constructors
    .def(py::init<>())
    .def(py::init<const GCP&>(), "rhs"_a)
    .def(py::init<const GCI&, const GCI&>(), "point"_a, "normal"_a)

    // Methods
    .def("minimumDistance", &GCP::minimumDistance)
    .def("coplanar", &GCP::coplanar)
    .def("parallel", &GCP::parallel)
    .def("valid", &GCP::valid)

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    .def(py::self >  GCI())
    .def(py::self <  GCI())
    .def(py::self >= GCI())
    .def(py::self <= GCI())

    // Attributes
    .def_property("point", &GCP::point, &GCP::setPoint)
    .def_property("normal", &GCP::normal, &GCP::setNormal)
    ;
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralNeighbor) {
  py::module m("SpheralNeighbor", "Spheral Neighbor module.");

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif
// #ifdef SPHERAL2D
//   dimensionBindings<Spheral::Dim<2>>(m, "2d");
// #endif
// #ifdef SPHERAL3D
//   dimensionBindings<Spheral::Dim<3>>(m, "3d");
// #endif

  return m.ptr();
}
