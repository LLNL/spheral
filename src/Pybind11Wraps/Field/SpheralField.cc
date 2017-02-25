#include <vector>
#include <string>

#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include "Geometry/Dimension.hh"
#include "Field/FieldBase.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "Utilities/FieldDataTypeTraits.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral::FieldSpace;

namespace {  // anonymous

//------------------------------------------------------------------------------
// Bind methods to Field objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename FieldObj, typename PB11Obj>
void fieldBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  using Spheral::NodeSpace::NodeList;

  obj

    // Constructors
    .def(py::init<std::string>(), "name"_a)
    .def(py::init<std::string, const FieldObj&>(), "name"_a, "field"_a)
    .def(py::init<std::string, const NodeList<Dimension>&>(), "name"_a, "nodeList"_a)
    .def(py::init<std::string, const NodeList<Dimension>&, typename FieldObj::value_type>(), "name"_a, "nodeList"_a, "value"_a)
    .def(py::init<const FieldObj&>(), "field"_a)

    // Methods
    .def("Zero", &FieldObj::Zero)
    .def("valid", &FieldObj::valid)
    .def("internalValues", &FieldObj::internalValues)
    .def("ghostValues", &FieldObj::ghostValues)
    .def("allValues", &FieldObj::allValues)

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    
    // Sequence operations
    .def("__getitem__", [](const FieldObj &s, size_t i) {
        if (i >= s.size()) throw py::index_error();
        return s[i];
      })
    .def("__setitem__", [](FieldObj &s, size_t i, typename FieldObj::value_type& v) {
        if (i >= s.size()) throw py::index_error();
        s[i] = v;
      })
    .def("__len__", [](const FieldObj& self) { return self.size(); })
    .def("__iter__", [](const FieldObj &s) { return py::make_iterator(s.begin(), s.end()); },
         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("__contains__", [](const FieldObj& self, const typename FieldObj::value_type& val) { for (const auto& elem: self) { if (elem == val) return true; } return false; })

    // We like in index with operator()
    .def("__call__", (typename FieldObj::value_type& (FieldObj::*)(int)) &FieldObj::operator(), py::is_operator(), py::return_value_policy::reference_internal)

    // FieldBase methods, since we currently don't expose the FieldBase to python.
    .def("size", &FieldObj::size)
  
    // Attributes
    .def_property_readonly("numElements", &FieldObj::numElements)
    .def_property_readonly("numInternalElements", &FieldObj::numInternalElements)

    ;
}

//------------------------------------------------------------------------------
// Bind numeric methods to Field objects, since only some Fields support such
// operations.
//------------------------------------------------------------------------------
template<typename Dimension, typename FieldObj, typename PB11Obj>
void fieldNumericBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  obj
    .def("applyMin", &FieldObj::applyMin)
    .def("applyMax", &FieldObj::applyMax)
    .def("sumElements", &FieldObj::sumElements)
    .def("min", &FieldObj::min)
    .def("max", &FieldObj::max)
    .def("localSumElements", &FieldObj::localSumElements)
    .def("localMin", &FieldObj::localMin)
    .def("localMax", &FieldObj::localMax)
    .def("string", (std::string (FieldObj::*)(const int) const) &FieldObj::string, "precision"_a=20)
    .def("string", (void (FieldObj::*)(const std::string&)) &FieldObj::string)
    
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self + typename FieldObj::value_type())
    .def(py::self - typename FieldObj::value_type())
    .def(py::self += typename FieldObj::value_type())
    .def(py::self -= typename FieldObj::value_type())

    .def(py::self *= float())
    .def(py::self /= float())
    .def(py::self *= double())
    .def(py::self /= double())

    .def(py::self <  py::self)
    .def(py::self >  py::self)
    .def(py::self <= py::self)
    .def(py::self >= py::self)

    ;
}

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  using namespace Spheral::FieldSpace;

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  //............................................................................
  // numeric Fields
  py::class_<Field<Dimension, int>> intfieldPB11(m, ("IntField" + suffix).c_str());
  py::class_<Field<Dimension, uint64_t>> ullfieldPB11(m, ("ULLField" + suffix).c_str());
  py::class_<Field<Dimension, Scalar>> scalarfieldPB11(m, ("ScalarField" + suffix).c_str());
  py::class_<Field<Dimension, Vector>> vectorfieldPB11(m, ("VectorField" + suffix).c_str());
  py::class_<Field<Dimension, Tensor>> tensorfieldPB11(m, ("TensorField" + suffix).c_str());
  py::class_<Field<Dimension, SymTensor>> symtensorfieldPB11(m, ("SymTensorField" + suffix).c_str());
  py::class_<Field<Dimension, ThirdRankTensor>> thirdranktensorfieldPB11(m, ("ThirdRankTensorField" + suffix).c_str());
  py::class_<Field<Dimension, FourthRankTensor>> fourthranktensorfieldPB11(m, ("FourthRankTensorField" + suffix).c_str());
  py::class_<Field<Dimension, FifthRankTensor>> fifthranktensorfieldPB11(m, ("FifthRankTensorField" + suffix).c_str());

  fieldBindings<Dimension, Field<Dimension, int>>(m, suffix, intfieldPB11);
  fieldBindings<Dimension, Field<Dimension, uint64_t>>(m, suffix, ullfieldPB11);
  fieldBindings<Dimension, Field<Dimension, Scalar>>(m, suffix, scalarfieldPB11);
  fieldBindings<Dimension, Field<Dimension, Vector>>(m, suffix, vectorfieldPB11);
  fieldBindings<Dimension, Field<Dimension, Tensor>>(m, suffix, tensorfieldPB11);
  fieldBindings<Dimension, Field<Dimension, SymTensor>>(m, suffix, symtensorfieldPB11);
  fieldBindings<Dimension, Field<Dimension, ThirdRankTensor>>(m, suffix, thirdranktensorfieldPB11);
  fieldBindings<Dimension, Field<Dimension, FourthRankTensor>>(m, suffix, fourthranktensorfieldPB11);
  fieldBindings<Dimension, Field<Dimension, FifthRankTensor>>(m, suffix, fifthranktensorfieldPB11);

  fieldNumericBindings<Dimension, Field<Dimension, int>>(m, suffix, intfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, uint64_t>>(m, suffix, ullfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, Scalar>>(m, suffix, scalarfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, Vector>>(m, suffix, vectorfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, Tensor>>(m, suffix, tensorfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, SymTensor>>(m, suffix, symtensorfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, ThirdRankTensor>>(m, suffix, thirdranktensorfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, FourthRankTensor>>(m, suffix, fourthranktensorfieldPB11);
  fieldNumericBindings<Dimension, Field<Dimension, FifthRankTensor>>(m, suffix, fifthranktensorfieldPB11);

  //............................................................................
  // non-numeric Fields
  py::class_<Field<Dimension, FacetedVolume>> facetedvolumefieldPB11(m, ("FacetedVolumeField" + suffix).c_str());
  py::class_<Field<Dimension, std::vector<Scalar>>> vectorscalarfieldPB11(m, ("VectorScalarField" + suffix).c_str());
  py::class_<Field<Dimension, std::vector<Vector>>> vectorvectorfieldPB11(m, ("VectorVectorField" + suffix).c_str());
  py::class_<Field<Dimension, std::vector<Tensor>>> vectortensorfieldPB11(m, ("VectorTensorField" + suffix).c_str());
  py::class_<Field<Dimension, std::vector<SymTensor>>> vectorsymtensorfieldPB11(m, ("VectorSymTensorField" + suffix).c_str());

  fieldBindings<Dimension, Field<Dimension, FacetedVolume>>(m, suffix, facetedvolumefieldPB11);
  fieldBindings<Dimension, Field<Dimension, std::vector<Scalar>>>(m, suffix, vectorscalarfieldPB11);
  fieldBindings<Dimension, Field<Dimension, std::vector<Vector>>>(m, suffix, vectorvectorfieldPB11);
  fieldBindings<Dimension, Field<Dimension, std::vector<Tensor>>>(m, suffix, vectortensorfieldPB11);
  fieldBindings<Dimension, Field<Dimension, std::vector<SymTensor>>>(m, suffix, vectorsymtensorfieldPB11);

  //............................................................................
  // STL containers
  // py::bind_vector<std::vector<NodeList<Dimension>*>>(m, "vector_of_NodeList" + suffix);
  // py::bind_vector<std::vector<FluidNodeList<Dimension>*>>(m, "vector_of_FluidNodeList" + suffix);

  // //............................................................................
  // // Methods
  // m.def("generateVoidNodes", &generateVoidNodes<Dimension>);
  // m.def("zerothNodalMoment", &nthNodalMoment<Dimension, 0U>);
  // m.def("firstNodalMoment", &nthNodalMoment<Dimension, 1>);
  // m.def("secondNodalMoment", &nthNodalMoment<Dimension, 2U>);
  // m.def("zerothAndFirstNodalMoments", &zerothAndFirstNodalMoments<Dimension>);

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralField) {
  py::module m("SpheralField", "Spheral Field module.");

  //............................................................................
  // FieldStorageType
  py::enum_<Spheral::FieldSpace::FieldStorageType>(m, "FieldStorageType")
    .value("Reference", Spheral::FieldSpace::Reference)
    .value("Copyo", Spheral::FieldSpace::Copy)
    .export_values();

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
