// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "Field/FieldBase.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Field/FieldListSet.hh"
#include "Utilities/FieldDataTypeTraits.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::FieldSpace;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
typedef Field<Dim<1>, Dim<1>::Scalar> ScalarField1d;
typedef Field<Dim<1>, Dim<1>::Vector> VectorField1d;
typedef Field<Dim<1>, Dim<1>::Tensor> TensorField1d;
typedef Field<Dim<1>, Dim<1>::SymTensor> SymTensorField1d;
PYBIND11_MAKE_OPAQUE(std::vector<ScalarField1d>);
PYBIND11_MAKE_OPAQUE(std::vector<VectorField1d>);
PYBIND11_MAKE_OPAQUE(std::vector<TensorField1d>);
PYBIND11_MAKE_OPAQUE(std::vector<SymTensorField1d>);

PYBIND11_MAKE_OPAQUE(std::vector<ScalarField1d*>);
PYBIND11_MAKE_OPAQUE(std::vector<VectorField1d*>);
PYBIND11_MAKE_OPAQUE(std::vector<TensorField1d*>);
PYBIND11_MAKE_OPAQUE(std::vector<SymTensorField1d*>);

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
typedef Field<Dim<2>, Dim<2>::Scalar> ScalarField2d;
typedef Field<Dim<2>, Dim<2>::Vector> VectorField2d;
typedef Field<Dim<2>, Dim<2>::Tensor> TensorField2d;
typedef Field<Dim<2>, Dim<2>::SymTensor> SymTensorField2d;
PYBIND11_MAKE_OPAQUE(std::vector<ScalarField2d>);
PYBIND11_MAKE_OPAQUE(std::vector<VectorField2d>);
PYBIND11_MAKE_OPAQUE(std::vector<TensorField2d>);
PYBIND11_MAKE_OPAQUE(std::vector<SymTensorField2d>);

PYBIND11_MAKE_OPAQUE(std::vector<ScalarField2d*>);
PYBIND11_MAKE_OPAQUE(std::vector<VectorField2d*>);
PYBIND11_MAKE_OPAQUE(std::vector<TensorField2d*>);
PYBIND11_MAKE_OPAQUE(std::vector<SymTensorField2d*>);

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
typedef Field<Dim<3>, Dim<3>::Scalar> ScalarField3d;
typedef Field<Dim<3>, Dim<3>::Vector> VectorField3d;
typedef Field<Dim<3>, Dim<3>::Tensor> TensorField3d;
typedef Field<Dim<3>, Dim<3>::SymTensor> SymTensorField3d;
PYBIND11_MAKE_OPAQUE(std::vector<ScalarField3d>);
PYBIND11_MAKE_OPAQUE(std::vector<VectorField3d>);
PYBIND11_MAKE_OPAQUE(std::vector<TensorField3d>);
PYBIND11_MAKE_OPAQUE(std::vector<SymTensorField3d>);

PYBIND11_MAKE_OPAQUE(std::vector<ScalarField3d*>);
PYBIND11_MAKE_OPAQUE(std::vector<VectorField3d*>);
PYBIND11_MAKE_OPAQUE(std::vector<TensorField3d*>);
PYBIND11_MAKE_OPAQUE(std::vector<SymTensorField3d*>);

namespace {  // anonymous

//------------------------------------------------------------------------------
// Bind methods to Field objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void fieldBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  using Spheral::NodeSpace::NodeList;

  obj

    // Constructors
    .def(py::init<std::string>(), "name"_a)
    .def(py::init<std::string, const Obj&>(), "name"_a, "field"_a)
    .def(py::init<std::string, const NodeList<Dimension>&>(), "name"_a, "nodeList"_a)
    .def(py::init<std::string, const NodeList<Dimension>&, typename Obj::value_type>(), "name"_a, "nodeList"_a, "value"_a)
    .def(py::init<const Obj&>(), "field"_a)

    // Methods
    .def("Zero", &Obj::Zero)
    .def("valid", &Obj::valid)
    .def("internalValues", &Obj::internalValues)
    .def("ghostValues", &Obj::ghostValues)
    .def("allValues", &Obj::allValues)

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    
    // Sequence operations
    .def("__getitem__", [](const Obj &s, size_t i) {
        if (i >= s.size()) throw py::index_error();
        return s[i];
      })
    .def("__setitem__", [](Obj &s, size_t i, typename Obj::value_type& v) {
        if (i >= s.size()) throw py::index_error();
        s[i] = v;
      })
    .def("__len__", [](const Obj& self) { return self.size(); })
    .def("__iter__", [](const Obj &s) { return py::make_iterator(s.begin(), s.end()); },
         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("__contains__", [](const Obj& self, const typename Obj::value_type& val) { for (const auto& elem: self) { if (elem == val) return true; } return false; })

    // We like in index with operator()
    .def("__call__", (typename Obj::value_type& (Obj::*)(int)) &Obj::operator(), py::is_operator(), py::return_value_policy::reference_internal)

    // FieldBase methods, since we currently don't expose the FieldBase to python.
    .def("size", &Obj::size)
  
    // Attributes
    .def_property_readonly("numElements", &Obj::numElements)
    .def_property_readonly("numInternalElements", &Obj::numInternalElements)

    ;
}

//------------------------------------------------------------------------------
// Bind numeric methods to Field objects, since only some Fields support such
// operations.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void fieldNumericBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  obj
    .def("applyMin", &Obj::applyMin)
    .def("applyMax", &Obj::applyMax)
    .def("sumElements", &Obj::sumElements)
    .def("min", &Obj::min)
    .def("max", &Obj::max)
    .def("localSumElements", &Obj::localSumElements)
    .def("localMin", &Obj::localMin)
    .def("localMax", &Obj::localMax)
    .def("string", (std::string (Obj::*)(const int) const) &Obj::string, "precision"_a=20)
    .def("string", (void (Obj::*)(const std::string&)) &Obj::string)
    
    .def(py::self + py::self)
    .def(py::self - py::self)
    .def(py::self += py::self)
    .def(py::self -= py::self)
    .def(py::self + typename Obj::value_type())
    .def(py::self - typename Obj::value_type())
    .def(py::self += typename Obj::value_type())
    .def(py::self -= typename Obj::value_type())

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
// Bind methods to FieldList objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void fieldlistBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  using Spheral::NodeSpace::NodeList;

  obj

    // Constructors
    .def(py::init<>())
    .def(py::init<FieldStorageType>(), "aStorageType"_a)

    // Methods
    .def("copyFields", &Obj::copyFields)
    .def("haveField", &Obj::haveField)
    .def("haveNodeList", &Obj::haveNodeList)
    .def("assignFields", &Obj::assignFields)
    .def("appendField", &Obj::appendField)
    .def("deleteField", &Obj::deleteField)
    .def("appendNewField", &Obj::appendNewField)
    .def("setMasterNodeLists", (void (Obj::*)(const Vector&, const SymTensor&) const) &Obj::setMasterNodeLists)
    .def("setMasterNodeLists", (void (Obj::*)(const Vector&) const) &Obj::setMasterNodeLists)
    .def("setRefineNodeLists", (void (Obj::*)(const Vector&, const SymTensor&) const) &Obj::setRefineNodeLists)
    .def("setRefineNodeLists", (void (Obj::*)(const Vector&) const) &Obj::setRefineNodeLists)
    .def("Zero", &Obj::Zero)
    .def("nodeListPtrs", &Obj::nodeListPtrs)
    .def("fieldForNodeList", [](const Obj& self, const NodeList<Dimension>& x) { return *(self.fieldForNodeList(x)); })

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    
    // Sequence operations
    .def("__getitem__", [](const Obj &s, size_t i) {
        if (i >= s.size()) throw py::index_error();
        return s[i];
      })
    .def("__setitem__", [](Obj &s, size_t i, typename Obj::value_type v) {
        if (i >= s.size()) throw py::index_error();
        *(s.begin() + i) = v;
      })
    .def("__len__", [](const Obj& self) { return self.size(); })
    .def("__iter__", [](const Obj &s) { return py::make_iterator(s.begin(), s.end()); },
         py::keep_alive<0, 1>() /* Essential: keep object alive while iterator exists */)
    .def("__contains__", [](const Obj& self, const typename Obj::value_type& val) { for (const auto& elem: self) { if (elem == val) return true; } return false; })

    // We like in index with operator()
    .def("__call__", (typename Obj::FieldDataType& (Obj::*)(unsigned, unsigned)) &Obj::operator(), py::is_operator(), py::return_value_policy::reference_internal)

    // FieldBase methods, since we currently don't expose the FieldBase to python.
    .def("size", &Obj::size)
  
    // Attributes
    .def_property_readonly("storageType", &Obj::storageType)
    .def_property_readonly("numFields", &Obj::numFields)
    .def_property_readonly("numNodes", &Obj::numNodes)
    .def_property_readonly("numInternalNodes", &Obj::numInternalNodes)
    .def_property_readonly("numGhostNodes", &Obj::numGhostNodes)

    ;
}

//------------------------------------------------------------------------------
// Bind numeric methods to FieldList objects, since only some Fields support such
// operations.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void fieldlistNumericBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  obj
    .def("applyMin", &Obj::applyMin)
    .def("applyMax", &Obj::applyMax)
    .def("sumElements", &Obj::sumElements)
    .def("min", &Obj::min)
    .def("max", &Obj::max)
    .def("localSumElements", &Obj::localSumElements)
    .def("localMin", &Obj::localMin)
    .def("localMax", &Obj::localMax)
    
    .def(py::self +  py::self)
    .def(py::self -  py::self)
    .def(py::self += py::self)
    .def(py::self -= py::self)

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

  // A couple of extra methods for SymTensor Fields.
  symtensorfieldPB11
    .def("applyScalarMin", &Field<Dimension, SymTensor>::applyScalarMin)
    .def("applyScalarMax", &Field<Dimension, SymTensor>::applyScalarMax)
    ;

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
  // numeric FieldLists
  py::class_<FieldList<Dimension, int>> intfieldlistPB11(m, ("IntFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, uint64_t>> ullfieldlistPB11(m, ("ULLFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, Scalar>> scalarfieldlistPB11(m, ("ScalarFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, Vector>> vectorfieldlistPB11(m, ("VectorFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, Tensor>> tensorfieldlistPB11(m, ("TensorFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, SymTensor>> symtensorfieldlistPB11(m, ("SymTensorFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, ThirdRankTensor>> thirdranktensorfieldlistPB11(m, ("ThirdRankTensorFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, FourthRankTensor>> fourthranktensorfieldlistPB11(m, ("FourthRankTensorFieldList" + suffix).c_str());
  py::class_<FieldList<Dimension, FifthRankTensor>> fifthranktensorfieldlistPB11(m, ("FifthRankTensorFieldList" + suffix).c_str());

  fieldlistBindings<Dimension, FieldList<Dimension, int>>(m, suffix, intfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, uint64_t>>(m, suffix, ullfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, Scalar>>(m, suffix, scalarfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, Vector>>(m, suffix, vectorfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, Tensor>>(m, suffix, tensorfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, SymTensor>>(m, suffix, symtensorfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, ThirdRankTensor>>(m, suffix, thirdranktensorfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, FourthRankTensor>>(m, suffix, fourthranktensorfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, FifthRankTensor>>(m, suffix, fifthranktensorfieldlistPB11);

  fieldlistNumericBindings<Dimension, FieldList<Dimension, int>>(m, suffix, intfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, uint64_t>>(m, suffix, ullfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, Scalar>>(m, suffix, scalarfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, Vector>>(m, suffix, vectorfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, Tensor>>(m, suffix, tensorfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, SymTensor>>(m, suffix, symtensorfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, ThirdRankTensor>>(m, suffix, thirdranktensorfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, FourthRankTensor>>(m, suffix, fourthranktensorfieldlistPB11);
  fieldlistNumericBindings<Dimension, FieldList<Dimension, FifthRankTensor>>(m, suffix, fifthranktensorfieldlistPB11);

  // A couple of extra methods for SymTensor FieldLists.
  symtensorfieldlistPB11
    .def("applyScalarMin", &FieldList<Dimension, SymTensor>::applyScalarMin)
    .def("applyScalarMax", &FieldList<Dimension, SymTensor>::applyScalarMax)
    ;

  //............................................................................
  // non-numeric Fields
  py::class_<FieldList<Dimension, FacetedVolume>> facetedvolumefieldlistPB11(m, ("FacetedVolumeField" + suffix).c_str());
  py::class_<FieldList<Dimension, std::vector<Scalar>>> vectorscalarfieldlistPB11(m, ("VectorScalarField" + suffix).c_str());
  py::class_<FieldList<Dimension, std::vector<Vector>>> vectorvectorfieldlistPB11(m, ("VectorVectorField" + suffix).c_str());
  py::class_<FieldList<Dimension, std::vector<Tensor>>> vectortensorfieldlistPB11(m, ("VectorTensorField" + suffix).c_str());
  py::class_<FieldList<Dimension, std::vector<SymTensor>>> vectorsymtensorfieldlistPB11(m, ("VectorSymTensorField" + suffix).c_str());

  fieldlistBindings<Dimension, FieldList<Dimension, FacetedVolume>>(m, suffix, facetedvolumefieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, std::vector<Scalar>>>(m, suffix, vectorscalarfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, std::vector<Vector>>>(m, suffix, vectorvectorfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, std::vector<Tensor>>>(m, suffix, vectortensorfieldlistPB11);
  fieldlistBindings<Dimension, FieldList<Dimension, std::vector<SymTensor>>>(m, suffix, vectorsymtensorfieldlistPB11);

  //............................................................................
  // FieldListSet
  py::class_<FieldListSet<Dimension>>(m, ("FieldListSet" + suffix).c_str())
    
    // Constructors
    .def(py::init<>())

    // Attributes
    .def_readwrite("ScalarFieldLists", &FieldListSet<Dimension>::ScalarFieldLists)
    .def_readwrite("VectorFieldLists", &FieldListSet<Dimension>::VectorFieldLists)
    .def_readwrite("TensorFieldLists", &FieldListSet<Dimension>::TensorFieldLists)
    .def_readwrite("SymTensorFieldLists", &FieldListSet<Dimension>::SymTensorFieldLists)
    ;

  //............................................................................
  // STL containers
  py::bind_vector<std::vector<Field<Dimension, Scalar>>>(m, "vector_of_ScalarField" + suffix);
  py::bind_vector<std::vector<Field<Dimension, Vector>>>(m, "vector_of_VectorField" + suffix);
  py::bind_vector<std::vector<Field<Dimension, Tensor>>>(m, "vector_of_TensorField" + suffix);
  py::bind_vector<std::vector<Field<Dimension, SymTensor>>>(m, "vector_of_SymTensorField" + suffix);
  py::bind_vector<std::vector<Field<Dimension, ThirdRankTensor>>>(m, "vector_of_ThirdRankTensorField" + suffix);
  py::bind_vector<std::vector<Field<Dimension, FourthRankTensor>>>(m, "vector_of_FourthRankTensorField" + suffix);
  py::bind_vector<std::vector<Field<Dimension, FifthRankTensor>>>(m, "vector_of_FifthRankTensorField" + suffix);

  py::bind_vector<std::vector<Field<Dimension, Scalar>*>>(m, "vector_of_ScalarFieldPtr" + suffix);
  py::bind_vector<std::vector<Field<Dimension, Vector>*>>(m, "vector_of_VectorFieldPtr" + suffix);
  py::bind_vector<std::vector<Field<Dimension, Tensor>*>>(m, "vector_of_TensorFieldPtr" + suffix);
  py::bind_vector<std::vector<Field<Dimension, SymTensor>*>>(m, "vector_of_SymTensorFieldPtr" + suffix);
  py::bind_vector<std::vector<Field<Dimension, ThirdRankTensor>*>>(m, "vector_of_ThirdRankTensorFieldPtr" + suffix);
  py::bind_vector<std::vector<Field<Dimension, FourthRankTensor>*>>(m, "vector_of_FourthRankTensorFieldPtr" + suffix);
  py::bind_vector<std::vector<Field<Dimension, FifthRankTensor>*>>(m, "vector_of_FifthRankTensorFieldPtr" + suffix);
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
    .value("Reference", Spheral::FieldSpace::FieldStorageType::ReferenceFields)
    .value("Copyo", Spheral::FieldSpace::FieldStorageType::CopyFields)
    .export_values();

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif
#ifdef SPHERAL2D
  dimensionBindings<Spheral::Dim<2>>(m, "2d");
#endif
#ifdef SPHERAL3D
  dimensionBindings<Spheral::Dim<3>>(m, "3d");
#endif

  return m.ptr();
}
