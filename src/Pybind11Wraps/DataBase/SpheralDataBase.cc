// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/StateBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Neighbor/ConnectivityMap.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::DataBaseSpace;

namespace {  // anonymous

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  using Spheral::NodeSpace::NodeList;
  using Spheral::FieldSpace::Field;
  using Spheral::FieldSpace::FieldList;
  using Spheral::FieldSpace::FieldBase;
  using Spheral::FieldSpace::FieldListBase;

  //............................................................................
  // StateBase
  typedef StateBase<Dimension> SB;
  py::class_<SB>(m, ("StateBase" + suffix).c_str())

    // Constructors
    .def(py::init<>())

    // Methods
    .def("registered", (bool (SB::*)(const typename SB::KeyType& key) const) &SB::registered, "key"_a)
    .def("registered", (bool (SB::*)(const FieldBase<Dimension>& field) const) &SB::registered, "field"_a)
    .def("registered", (bool (SB::*)(const FieldListBase<Dimension>& fieldList) const) &SB::registered, "fieldList"_a)
    .def("fieldNameRegistered", &SB::fieldNameRegistered)
    .def("enroll", (void (SB::*)(FieldBase<Dimension>&)) &SB::enroll, "field"_a)
    .def("enroll", (void (SB::*)(FieldListBase<Dimension>&)) &SB::enroll, "fieldList"_a)
    .def("enrollMesh", &SB::enrollMesh)
    .def("enroll", (void (SB::*)(const typename SB::FieldName&, std::vector<typename Dimension::Scalar>&)) &SB::enroll, "key"_a, "vec"_a)
    .def("intField", (Field<Dimension, int>& (SB::*)(const typename SB::KeyType&, const int&) const) &SB::field, "key"_a, "dummy"_a)
    .def("scalarField", (Field<Dimension, double>& (SB::*)(const typename SB::KeyType&, const double&) const) &SB::field, "key"_a, "dummy"_a)
    .def("vectorField", (Field<Dimension, Vector>& (SB::*)(const typename SB::KeyType&, const Vector&) const) &SB::field, "key"_a, "dummy"_a)
    .def("tensorField", (Field<Dimension, Tensor>& (SB::*)(const typename SB::KeyType&, const Tensor&) const) &SB::field, "key"_a, "dummy"_a)
    .def("symTensorField", (Field<Dimension, SymTensor>& (SB::*)(const typename SB::KeyType&, const SymTensor&) const) &SB::field, "key"_a, "dummy"_a)
    .def("intFields", (FieldList<Dimension, int> (SB::*)(const typename std::string&, const int&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("scalarFields", (FieldList<Dimension, double> (SB::*)(const typename std::string&, const double&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("vectorFields", (FieldList<Dimension, Vector> (SB::*)(const typename std::string&, const Vector&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("tensorFields", (FieldList<Dimension, Tensor> (SB::*)(const typename std::string&, const Tensor&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("symTensorFields", (FieldList<Dimension, SymTensor> (SB::*)(const typename std::string&, const SymTensor&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("thirdRankTensorFields", (FieldList<Dimension, ThirdRankTensor> (SB::*)(const typename std::string&, const ThirdRankTensor&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("vectorDoubleFields", (FieldList<Dimension, std::vector<double>> (SB::*)(const typename std::string&, const std::vector<double>&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("vectorVectorFields", (FieldList<Dimension, std::vector<Vector>> (SB::*)(const typename std::string&, const std::vector<Vector>&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("vectorSymTensorFields", (FieldList<Dimension, std::vector<SymTensor>> (SB::*)(const typename std::string&, const std::vector<SymTensor>&) const) &SB::fields, "name"_a, "dummy"_a)

    // Comparisons
    .def(py::self == py::self)
    ;

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralDataBase) {
  py::module m("SpheralDataBase", "Spheral DataBase module.");

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
