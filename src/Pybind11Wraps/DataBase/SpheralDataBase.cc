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
  typedef typename Dimension::FourthRankTensor FourthRankTensor;
  typedef typename Dimension::FifthRankTensor FifthRankTensor;
  typedef typename Dimension::FacetedVolume FacetedVolume;

  using Spheral::NodeSpace::NodeList;
  using Spheral::FieldSpace::Field;
  using Spheral::FieldSpace::FieldList;
  using Spheral::FieldSpace::FieldBase;
  using Spheral::FieldSpace::FieldListBase;

  typedef StateBase<Dimension> SB;
  typedef State<Dimension> ST;
  typedef StateDerivatives<Dimension> SD;
  typedef DataBase<Dimension> DB;

  //............................................................................
  // StateBase
  py::class_<SB>(m, ("StateBase" + suffix).c_str())

    // Constructors
    .def(py::init<>())
    .def(py::init<const SB&>())

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
    .def("allIntFields", (std::vector<Field<Dimension, int>*> (SB::*)(const int&) const) &SB::allFields, "dummy"_a=0)
    .def("allScalarFields", (std::vector<Field<Dimension, double>*> (SB::*)(const double&) const) &SB::allFields, "dummy"_a=0.0)
    .def("allVectorFields", (std::vector<Field<Dimension, Vector>*> (SB::*)(const Vector&) const) &SB::allFields, "dummy"_a=Vector::zero)
    .def("allTensorFields", (std::vector<Field<Dimension, Tensor>*> (SB::*)(const Tensor&) const) &SB::allFields, "dummy"_a=Tensor::zero)
    .def("allSymTensorFields", (std::vector<Field<Dimension, SymTensor>*> (SB::*)(const SymTensor&) const) &SB::allFields, "dummy"_a=SymTensor::zero)

    // These are straight "field", "fields", and "allFields" overloaded versions of above, more like we do in C++.  Maybe with pybind11 this will work better?
    .def("field", (Field<Dimension, int>& (SB::*)(const typename SB::KeyType&, const int&) const) &SB::field, "key"_a, "dummy"_a)
    .def("field", (Field<Dimension, double>& (SB::*)(const typename SB::KeyType&, const double&) const) &SB::field, "key"_a, "dummy"_a)
    .def("field", (Field<Dimension, Vector>& (SB::*)(const typename SB::KeyType&, const Vector&) const) &SB::field, "key"_a, "dummy"_a)
    .def("field", (Field<Dimension, Tensor>& (SB::*)(const typename SB::KeyType&, const Tensor&) const) &SB::field, "key"_a, "dummy"_a)
    .def("field", (Field<Dimension, SymTensor>& (SB::*)(const typename SB::KeyType&, const SymTensor&) const) &SB::field, "key"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, int> (SB::*)(const typename std::string&, const int&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, double> (SB::*)(const typename std::string&, const double&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, Vector> (SB::*)(const typename std::string&, const Vector&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, Tensor> (SB::*)(const typename std::string&, const Tensor&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, SymTensor> (SB::*)(const typename std::string&, const SymTensor&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, ThirdRankTensor> (SB::*)(const typename std::string&, const ThirdRankTensor&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, std::vector<double>> (SB::*)(const typename std::string&, const std::vector<double>&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, std::vector<Vector>> (SB::*)(const typename std::string&, const std::vector<Vector>&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("fields", (FieldList<Dimension, std::vector<SymTensor>> (SB::*)(const typename std::string&, const std::vector<SymTensor>&) const) &SB::fields, "name"_a, "dummy"_a)
    .def("allFields", (std::vector<Field<Dimension, int>*> (SB::*)(const int&) const) &SB::allFields, "dummy"_a)
    .def("allFields", (std::vector<Field<Dimension, double>*> (SB::*)(const double&) const) &SB::allFields, "dummy"_a)
    .def("allFields", (std::vector<Field<Dimension, Vector>*> (SB::*)(const Vector&) const) &SB::allFields, "dummy"_a)
    .def("allFields", (std::vector<Field<Dimension, Tensor>*> (SB::*)(const Tensor&) const) &SB::allFields, "dummy"_a)
    .def("allFields", (std::vector<Field<Dimension, SymTensor>*> (SB::*)(const SymTensor&) const) &SB::allFields, "dummy"_a)

    .def("keys", &SB::keys)
    .def("fieldKeys", &SB::fieldKeys)
    .def("meshRegistered", &SB::meshRegistered)
    .def("mesh", (typename SB::MeshType& (SB::*)()) &SB::mesh)
    .def("assign", &SB::assign)
    .def("copyState", &SB::copyState)

    // Static methods
    .def_static("key", (typename SB::KeyType (*)(const FieldBase<Dimension>&)) &SB::key)
    .def_static("key", (typename SB::KeyType (*)(const FieldListBase<Dimension>&)) &SB::key)
    .def_static("buildFieldKey", &SB::buildFieldKey)
    .def_static("splitFieldKey", &SB::splitFieldKey)

    // Comparisons
    .def(py::self == py::self)
    ;

  //............................................................................
  // State
  py::class_<ST, SB> pyST(m, ("State" + suffix).c_str());

  // Constructors
  pyST.def(py::init<>());
  pyST.def(py::init<DB&, typename ST::PackageList&>(), "dataBase"_a, "packages"_a);
  pyST.def(py::init<const ST&>());

  // Methods
  pyST.def("update", &ST::update);
  pyST.def("policyKeys", &ST::policyKeys);

  // Comparisons
  pyST.def(py::self == py::self);

  //............................................................................
  // StateDerivatives
  py::class_<SD, SB> pySD(m, ("StateDerivatives" + suffix).c_str());

  // Constructors
  pySD.def(py::init<>());
  pySD.def(py::init<DB&, typename SD::PackageList&>(), "dataBase"_a, "packages"_a);
  pySD.def(py::init<const SD&>());

  // Methods
  pySD.def("initializeNodePairInformation", &SD::initializeNodePairInformation);
  pySD.def("calculatedNodePairsSymmetric", &SD::calculatedNodePairsSymmetric);
  pySD.def("Zero", &SD::Zero);

  // Comparisons
  pySD.def(py::self == py::self);

  //............................................................................
  // DataBase
  py::class_<DB> pyDB(m, ("DataBase" + suffix).c_str());

  // Constructors
  pyDB.def(py::init<>());

  // Methods
  pyDB.def("reinitializeNeighbors", &DB::reinitializeNeighbors);
  pyDB.def("updateConnectivityMap", &DB::updateConnectivityMap);
  pyDB.def("patchConnectivityMap", &DB::patchConnectivityMap);
  pyDB.def("connectivityMap", (const NeighborSpace::ConnectivityMap<Dimension>& (DB::*)(const bool) const) &DB::connectivityMap, "computeGhostConnectivity"_a=false);
  pyDB.def("appendNodeList", (void (DB::*)(NodeSpace::SolidNodeList<Dimension>&)) &DB::appendNodeList, "nodeList"_a);
  pyDB.def("appendNodeList", (void (DB::*)(NodeSpace::FluidNodeList<Dimension>&)) &DB::appendNodeList, "nodeList"_a);
  pyDB.def("appendNodeList", (void (DB::*)(NodeSpace::NodeList<Dimension>&)) &DB::appendNodeList, "nodeList"_a);
  pyDB.def("deleteNodeList", (void (DB::*)(NodeSpace::SolidNodeList<Dimension>&)) &DB::deleteNodeList, "nodeList"_a);
  pyDB.def("deleteNodeList", (void (DB::*)(NodeSpace::FluidNodeList<Dimension>&)) &DB::deleteNodeList, "nodeList"_a);
  pyDB.def("deleteNodeList", (void (DB::*)(NodeSpace::NodeList<Dimension>&)) &DB::deleteNodeList, "nodeList"_a);
  pyDB.def("haveNodeList", &DB::haveNodeList);
  pyDB.def("nodeListPtrs", &DB::nodeListPtrs);
  pyDB.def("fluidNodeListPtrs", &DB::fluidNodeListPtrs);
  pyDB.def("solidNodeListPtrs", &DB::solidNodeListPtrs);
  pyDB.def("setMasterNodeLists", &DB::setMasterNodeLists, "position"_a, "H"_a, "masterLists"_a, "coarseNeighbors"_a);
  pyDB.def("setMasterFluidNodeLists", &DB::setMasterFluidNodeLists, "position"_a, "H"_a, "masterLists"_a, "coarseNeighbors"_a);
  pyDB.def("setRefineNodeLists", &DB::setRefineNodeLists, "position"_a, "H"_a, "coarseNeighbors"_a, "refineNeighbors"_a);
  pyDB.def("setRefineFluidNodeLists", &DB::setRefineFluidNodeLists, "position"_a, "H"_a, "coarseNeighbors"_a, "refineNeighbors"_a);
  pyDB.def("globalMass", &DB::globalMass);
  pyDB.def("globalPosition", &DB::globalPosition);
  pyDB.def("globalVelocity", &DB::globalVelocity);
  pyDB.def("globalHfield", &DB::globalHfield);
  pyDB.def("globalWork", &DB::globalWork);
  pyDB.def("fluidMass", &DB::fluidMass);
  pyDB.def("fluidPosition", &DB::fluidPosition);
  pyDB.def("fluidVelocity", &DB::fluidVelocity);
  pyDB.def("fluidMassDensity", &DB::fluidMassDensity);
  pyDB.def("fluidSpecificThermalEnergy", &DB::fluidSpecificThermalEnergy);
  pyDB.def("fluidHfield", &DB::fluidHfield);
  pyDB.def("fluidWork", &DB::fluidWork);
  pyDB.def("globalNodeExtent", &DB::globalNodeExtent);
  pyDB.def("fluidNodeExtent", &DB::fluidNodeExtent);
  pyDB.def("globalHinverse", &DB::globalHinverse);
  pyDB.def("fluidHinverse", &DB::fluidHinverse);
  pyDB.def("fluidPressure", &DB::fluidPressure);
  pyDB.def("fluidTemperature", &DB::fluidTemperature);
  pyDB.def("fluidSoundSpeed", &DB::fluidSoundSpeed);
  pyDB.def("fluidVolume", &DB::fluidVolume);
  pyDB.def("fluidGamma", &DB::fluidGamma);
  pyDB.def("fluidEntropy", &DB::fluidEntropy);
  pyDB.def("fluidLinearMomentum", &DB::fluidLinearMomentum);
  pyDB.def("fluidTotalEnergy", &DB::fluidTotalEnergy);

  pyDB.def("newGlobalIntFieldList", (FieldList<Dimension, int> (DB::*)(const int, const std::string) const) &DB::newGlobalFieldList, "value"_a=0, "name"_a="unnamed field list");
  pyDB.def("newGlobalScalarFieldList", (FieldList<Dimension, double> (DB::*)(const double, const std::string) const) &DB::newGlobalFieldList, "value"_a=0.0, "name"_a="unnamed field list");
  pyDB.def("newGlobalVectorFieldList", (FieldList<Dimension, Vector> (DB::*)(const Vector, const std::string) const) &DB::newGlobalFieldList, "value"_a=Vector::zero, "name"_a="unnamed field list");
  pyDB.def("newGlobalTensorFieldList", (FieldList<Dimension, Tensor> (DB::*)(const Tensor, const std::string) const) &DB::newGlobalFieldList, "value"_a=Tensor::zero, "name"_a="unnamed field list");
  pyDB.def("newGlobalSymTensorFieldList", (FieldList<Dimension, SymTensor> (DB::*)(const SymTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a=SymTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newGlobalThirdRankTensorFieldList", (FieldList<Dimension, ThirdRankTensor> (DB::*)(const ThirdRankTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a=ThirdRankTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newGlobalFourthRankTensorFieldList", (FieldList<Dimension, FourthRankTensor> (DB::*)(const FourthRankTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a=FourthRankTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newGlobalFifthRankTensorFieldList", (FieldList<Dimension, FifthRankTensor> (DB::*)(const FifthRankTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a=FifthRankTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newGlobalFacetedVolumeFieldList", (FieldList<Dimension, FacetedVolume> (DB::*)(const FacetedVolume, const std::string) const) &DB::newGlobalFieldList, "value"_a=FacetedVolume(), "name"_a="unnamed field list");
  pyDB.def("newGlobalVectorOfDoubleFieldList", (FieldList<Dimension, std::vector<double>> (DB::*)(const std::vector<double>, const std::string) const) &DB::newGlobalFieldList, "value"_a=std::vector<double>(), "name"_a="unnamed field list");

  pyDB.def("newFluidIntFieldList", (FieldList<Dimension, int> (DB::*)(const int, const std::string) const) &DB::newFluidFieldList, "value"_a=0, "name"_a="unnamed field list");
  pyDB.def("newFluidScalarFieldList", (FieldList<Dimension, double> (DB::*)(const double, const std::string) const) &DB::newFluidFieldList, "value"_a=0.0, "name"_a="unnamed field list");
  pyDB.def("newFluidVectorFieldList", (FieldList<Dimension, Vector> (DB::*)(const Vector, const std::string) const) &DB::newFluidFieldList, "value"_a=Vector::zero, "name"_a="unnamed field list");
  pyDB.def("newFluidTensorFieldList", (FieldList<Dimension, Tensor> (DB::*)(const Tensor, const std::string) const) &DB::newFluidFieldList, "value"_a=Tensor::zero, "name"_a="unnamed field list");
  pyDB.def("newFluidSymTensorFieldList", (FieldList<Dimension, SymTensor> (DB::*)(const SymTensor, const std::string) const) &DB::newFluidFieldList, "value"_a=SymTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newFluidThirdRankTensorFieldList", (FieldList<Dimension, ThirdRankTensor> (DB::*)(const ThirdRankTensor, const std::string) const) &DB::newFluidFieldList, "value"_a=ThirdRankTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newFluidFourthRankTensorFieldList", (FieldList<Dimension, FourthRankTensor> (DB::*)(const FourthRankTensor, const std::string) const) &DB::newFluidFieldList, "value"_a=FourthRankTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newFluidFifthRankTensorFieldList", (FieldList<Dimension, FifthRankTensor> (DB::*)(const FifthRankTensor, const std::string) const) &DB::newFluidFieldList, "value"_a=FifthRankTensor::zero, "name"_a="unnamed field list");
  pyDB.def("newFluidFacetedVolumeFieldList", (FieldList<Dimension, FacetedVolume> (DB::*)(const FacetedVolume, const std::string) const) &DB::newFluidFieldList, "value"_a=FacetedVolume(), "name"_a="unnamed field list");
  pyDB.def("newFluidVectorOfDoubleFieldList", (FieldList<Dimension, std::vector<double>> (DB::*)(const std::vector<double>, const std::string) const) &DB::newFluidFieldList, "value"_a=std::vector<double>(), "name"_a="unnamed field list");

  // Same as above without encoding the return type in method name
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, int> (DB::*)(const int, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, double> (DB::*)(const double, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, Vector> (DB::*)(const Vector, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, Tensor> (DB::*)(const Tensor, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, SymTensor> (DB::*)(const SymTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, ThirdRankTensor> (DB::*)(const ThirdRankTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, FourthRankTensor> (DB::*)(const FourthRankTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, FifthRankTensor> (DB::*)(const FifthRankTensor, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, FacetedVolume> (DB::*)(const FacetedVolume, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newGlobalFieldList", (FieldList<Dimension, std::vector<double>> (DB::*)(const std::vector<double>, const std::string) const) &DB::newGlobalFieldList, "value"_a, "name"_a="unnamed field list");

  pyDB.def("newFluidFieldList", (FieldList<Dimension, int> (DB::*)(const int, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, double> (DB::*)(const double, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, Vector> (DB::*)(const Vector, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, Tensor> (DB::*)(const Tensor, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, SymTensor> (DB::*)(const SymTensor, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, ThirdRankTensor> (DB::*)(const ThirdRankTensor, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, FourthRankTensor> (DB::*)(const FourthRankTensor, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, FifthRankTensor> (DB::*)(const FifthRankTensor, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, FacetedVolume> (DB::*)(const FacetedVolume, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");
  pyDB.def("newFluidFieldList", (FieldList<Dimension, std::vector<double>> (DB::*)(const std::vector<double>, const std::string) const) &DB::newFluidFieldList, "value"_a, "name"_a="unnamed field list");

  // Attributes
  pyDB.def_property_readonly("numNodeLists", &DB::numNodeLists);
  pyDB.def_property_readonly("numFluidNodeLists", &DB::numFluidNodeLists);
  pyDB.def_property_readonly("numSolidNodeLists", &DB::numSolidNodeLists);
  pyDB.def_property_readonly("numInternalNodes", &DB::numInternalNodes);
  pyDB.def_property_readonly("numGhostNodes", &DB::numGhostNodes);
  pyDB.def_property_readonly("numNodes", &DB::numNodes);
  pyDB.def_property_readonly("globalNumInternalNodes", &DB::globalNumInternalNodes);
  pyDB.def_property_readonly("globalNumGhostNodes", &DB::globalNumGhostNodes);
  pyDB.def_property_readonly("globalNumNodes", &DB::globalNumNodes);

  pyDB.def_property_readonly("globalMass", &DB::globalMass);
  pyDB.def_property_readonly("globalPosition", &DB::globalPosition);
  pyDB.def_property_readonly("globalVelocity", &DB::globalVelocity);
  pyDB.def_property_readonly("globalHfield", &DB::globalHfield);
  pyDB.def_property_readonly("globalWork", &DB::globalWork);

  pyDB.def_property_readonly("fluidMass", &DB::fluidMass);
  pyDB.def_property_readonly("fluidPosition", &DB::fluidPosition);
  pyDB.def_property_readonly("fluidVelocity", &DB::fluidVelocity);
  pyDB.def_property_readonly("fluidMassDensity", &DB::fluidMassDensity);
  pyDB.def_property_readonly("fluidSpecificThermalEnergy", &DB::fluidSpecificThermalEnergy);
  pyDB.def_property_readonly("fluidHfield", &DB::fluidHfield);
  pyDB.def_property_readonly("fluidWork", &DB::fluidWork);

  pyDB.def_property_readonly("globalNodeExtent", &DB::globalNodeExtent);
  pyDB.def_property_readonly("fluidNodeExtent", &DB::fluidNodeExtent);
  pyDB.def_property_readonly("numNeighbors", &DB::numNeighbors);
  pyDB.def_property_readonly("maxKernelExtent", &DB::maxKernelExtent);

  // For backwards compatibility, we provide attribute-like access to the static attributes.
  pyDB.def_property_readonly("nDim", []() { return DB::nDim; });
  pyDB.def_property_readonly("isRZ", []() { return DB::isRZ; });

  // Static attributes
  pyDB.def_readonly_static("nDim", &DB::nDim);
  pyDB.def_readonly_static("isRZ", &DB::isRZ);
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralDataBase, m) {
  m.doc() = "Spheral DataBase module.";

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
}
