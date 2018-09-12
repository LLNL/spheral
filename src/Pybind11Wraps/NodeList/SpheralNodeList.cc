// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "NodeList/NodeList.hh"
#include "NodeList/FluidNodeList.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "NodeList/FixedSmoothingScale.hh"
#include "NodeList/SPHSmoothingScale.hh"
#include "NodeList/ASPHSmoothingScale.hh"
#include "NodeList/generateVoidNodes.hh"
#include "NodeList/nthNodalMoment.hh"
#include "Kernel/TableKernel.hh"
#include "FileIO/FileIO.hh"
#include "Material/EquationOfState.hh"
#include "Neighbor/Neighbor.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Mesh/Mesh.hh"

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Spheral::NodeSpace;

namespace Spheral {
namespace NodeSpace {

//------------------------------------------------------------------------------
// PyNodeList
//------------------------------------------------------------------------------
template<class NodeListBase>
class PyNodeList: public NodeListBase {
public:
  using NodeListBase::NodeListBase;  // inherit constructors

  virtual void deleteNodes(const std::vector<int>& nodeIDs) override {
    PYBIND11_OVERLOAD(void,         // Return type
                      NodeListBase, // Parent class
                      deleteNodes,  // name of method
                      nodeIDs       // arguments
      );
  }

  virtual std::list<std::vector<char>> packNodeFieldValues(const std::vector<int>& nodeIDs) const override {
    PYBIND11_OVERLOAD(std::list<std::vector<char>>,  // Return type
                      NodeListBase,                  // Parent class
                      packNodeFieldValues,           // name of method
                      nodeIDs                        // arguments
      );
  }

  virtual void appendInternalNodes(const int numNewNodes, const std::list<std::vector<char>>& packedFieldValues) override {
    PYBIND11_OVERLOAD(void,                        // Return type
                      NodeListBase,                // Parent class
                      appendInternalNodes,         // name of method
                      numNewNodes,                 // arguments
                      packedFieldValues
      );
  }

  virtual void reorderNodes(const std::vector<int>& newOrdering) override {
    PYBIND11_OVERLOAD(void,                        // Return type
                      NodeListBase,                // Parent class
                      reorderNodes,                // name of method
                      newOrdering                  // arguments
      );
  }

  virtual std::string label() const override {
    PYBIND11_OVERLOAD(std::string,          // Return type
                      NodeListBase,         // Parent class
                      label                 // name of method
      );
  }

  virtual void dumpState(FileIOSpace::FileIO& file, const std::string& pathName) const override {
    PYBIND11_OVERLOAD(void,          // Return type
                      NodeListBase,  // Parent class
                      dumpState,     // name of method
                      file,          // arguments
                      pathName
      );
  }

  virtual void restoreState(const FileIOSpace::FileIO& file, const std::string& pathName) override {
    PYBIND11_OVERLOAD(void,                 // Return type
                      NodeListBase,         // Parent class
                      restoreState,         // name of method
                      file,                 // arguments
                      pathName
      );
  }

};

//------------------------------------------------------------------------------
// PyFluidNodeList
//------------------------------------------------------------------------------
template<typename Dimension, class FluidNodeListBase>
class PyFluidNodeList: public PyNodeList<FluidNodeListBase> {
public:
  using PyNodeList<FluidNodeListBase>::PyNodeList;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  virtual void pressure(FieldSpace::Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD(void,              // Return type
                      FluidNodeListBase, // Parent class
                      pressure,          // name of method
                      field              // arguments
      );
  }

  virtual void temperature(FieldSpace::Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD(void,              // Return type
                      FluidNodeListBase, // Parent class
                      temperature,       // name of method
                      field              // arguments
      );
  }

  virtual void soundSpeed(FieldSpace::Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD(void,              // Return type
                      FluidNodeListBase, // Parent class
                      soundSpeed,        // name of method
                      field              // arguments
      );
  }

  virtual void volume(FieldSpace::Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD(void,              // Return type
                      FluidNodeListBase, // Parent class
                      volume,            // name of method
                      field              // arguments
      );
  }

  virtual void linearMomentum(FieldSpace::Field<Dimension, Vector>& field) const override {
    PYBIND11_OVERLOAD(void,              // Return type
                      FluidNodeListBase, // Parent class
                      linearMomentum,    // name of method
                      field              // arguments
      );
  }

  virtual void totalEnergy(FieldSpace::Field<Dimension, Scalar>& field) const override {
    PYBIND11_OVERLOAD(void,              // Return type
                      FluidNodeListBase, // Parent class
                      totalEnergy,       // name of method
                      field              // arguments
      );
  }
};

//------------------------------------------------------------------------------
// PySmoothingScale
// Common methods for SmoothingScaleBase, SPHSmoothingScale, & ASPHSmoothingScale
//------------------------------------------------------------------------------
template<typename Dimension, class SSBase>
class PySmoothingScale: public SSBase {
public:
  using SSBase::SSBase;  // inherit constructors

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  virtual SymTensor smoothingScaleDerivative(const SymTensor& H,
                                             const Vector& pos,
                                             const Tensor& DvDx,
                                             const Scalar hmin,
                                             const Scalar hmax,
                                             const Scalar hminratio,
                                             const Scalar nPerh) const override {
    PYBIND11_OVERLOAD(SymTensor,                                 // Return type
                      SSBase,                                    // Parent class
                      smoothingScaleDerivative,                  // name of method
                      H, pos, DvDx, hmin, hmax, hminratio, nPerh // arguments
      );
  }

  virtual SymTensor newSmoothingScale(const SymTensor& H,
                                      const Vector& pos,
                                      const Scalar zerothMoment,
                                      const SymTensor& secondMoment,
                                      const KernelSpace::TableKernel<Dimension>& W,
                                      const Scalar hmin,
                                      const Scalar hmax,
                                      const Scalar hminratio,
                                      const Scalar nPerh,
                                      const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                      const unsigned nodeListi,
                                      const unsigned i) const override {
    PYBIND11_OVERLOAD(SymTensor,              // Return type
                      SSBase,                 // Parent class
                      newSmoothingScale,      // name of method
                      H, pos, zerothMoment, secondMoment, W, hmin, hmax, hminratio, nPerh, connectivityMap, nodeListi, i     // arguments
      );
  }

  virtual SymTensor idealSmoothingScale(const SymTensor& H,
                                        const Vector& pos,
                                        const Scalar zerothMoment,
                                        const SymTensor& secondMoment,
                                        const KernelSpace::TableKernel<Dimension>& W,
                                        const Scalar hmin,
                                        const Scalar hmax,
                                        const Scalar hminratio,
                                        const Scalar nPerh,
                                        const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                                        const unsigned nodeListi,
                                        const unsigned i) const override {
    PYBIND11_OVERLOAD(SymTensor,              // Return type
                      SSBase,                 // Parent class
                      idealSmoothingScale,    // name of method
                      H, pos, zerothMoment, secondMoment, W, hmin, hmax, hminratio, nPerh, connectivityMap, nodeListi, i     // arguments
      );
  }

  virtual SymTensor idealSmoothingScale(const SymTensor& H,
                                        const Spheral::MeshSpace::Mesh<Dimension>& mesh,
                                        const typename Spheral::MeshSpace::Mesh<Dimension>::Zone& zone,
                                        const Scalar hmin,
                                        const Scalar hmax,
                                        const Scalar hminratio,
                                        const Scalar nPerh) const override {
    PYBIND11_OVERLOAD(SymTensor,              // Return type
                      SSBase,                 // Parent class
                      idealSmoothingScale,    // name of method
                      H, mesh, zone, hmin, hmax, hminratio, nPerh     // arguments
      );
  }
};

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to nthNodalMoment
//------------------------------------------------------------------------------
template<typename Dimension, unsigned moment>
inline
FieldSpace::FieldList<Dimension, typename MomentTraits<Dimension, moment>::Moment>
nthNodalMoment(const std::vector<NodeList<Dimension>*>& nodeLists,
               const KernelSpace::TableKernel<Dimension>& W,
               const bool renormalize) {
  return nthNodalMoment<Dimension, typename std::vector<NodeList<Dimension>*>::const_iterator, moment>
    (nodeLists.begin(), nodeLists.end(), W, renormalize);
}

//------------------------------------------------------------------------------
// Provide a non-iterator based interface to zerothAndFirstNodalMoments
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
zerothAndFirstNodalMoments(const std::vector<NodeList<Dimension>*>& nodeLists,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const bool useGradientAsKernel,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& zerothMoment,
                           FieldSpace::FieldList<Dimension, typename Dimension::Vector>& firstMoment) {
  zerothAndFirstNodalMoments<Dimension, typename std::vector<NodeList<Dimension>*>::const_iterator>
    (nodeLists.begin(), nodeLists.end(), W, useGradientAsKernel, zerothMoment, firstMoment);
}

}
}

namespace {  // anonymous

//------------------------------------------------------------------------------
// Bind methods to Smoothing scale objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename SSObj, typename PB11Obj>
void smoothingScaleObjectBindings(py::module& m, const std::string suffix, PB11Obj& obj) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  obj

    // Constructors
    .def(py::init<>())

    // Methods
    .def("newSmoothingScaleAndDerivative", &SSObj::newSmoothingScaleAndDerivative)

    // Virtual methods
    .def("smoothingScaleDerivative", &SSObj::smoothingScaleDerivative)
    .def("newSmoothingScale", &SSObj::newSmoothingScale)
    .def("idealSmoothingScale",
         (SymTensor (SSObj::*)(const SymTensor&,
                               const Vector&,
                               const Scalar,
                               const SymTensor&,
                               const Spheral::KernelSpace::TableKernel<Dimension>&,
                               const Scalar,
                               const Scalar,
                               const Scalar,
                               const Scalar,
                               const Spheral::NeighborSpace::ConnectivityMap<Dimension>&,
                               const unsigned,
                               const unsigned) const) &SSObj::idealSmoothingScale)
    .def("idealSmoothingScale",
         (SymTensor (SSObj::*)(const SymTensor&,
                               const Spheral::MeshSpace::Mesh<Dimension>&,
                               const typename Spheral::MeshSpace::Mesh<Dimension>::Zone&,
                               const Scalar,
                               const Scalar,
                               const Scalar,
                               const Scalar) const) &SSObj::idealSmoothingScale)
    ;
}

//------------------------------------------------------------------------------
// Per dimension bindings.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  using namespace Spheral::FieldSpace;
  using Spheral::EquationOfState;

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;

  //............................................................................
  // NodeListRegistrar
  py::class_<Spheral::NodeListRegistrar<Dimension>,
             std::unique_ptr<Spheral::NodeListRegistrar<Dimension>, py::nodelete>>(m, ("NodeListRegistrar" + suffix).c_str(), py::metaclass())

    // Attributes
    .def_property_readonly_static("instance", &Spheral::NodeListRegistrar<Dimension>::instance)
    .def_property_readonly("numNodeLists", &Spheral::NodeListRegistrar<Dimension>::numNodeLists)
    .def_property_readonly("numFluidNodeLists", &Spheral::NodeListRegistrar<Dimension>::numFluidNodeLists)
    .def_property("domainDecompositionIndependent",
                  (bool (Spheral::NodeListRegistrar<Dimension>::*)() const) &Spheral::NodeListRegistrar<Dimension>::domainDecompositionIndependent,
                  (void (Spheral::NodeListRegistrar<Dimension>::*)(const bool)) &Spheral::NodeListRegistrar<Dimension>::domainDecompositionIndependent)
    .def_property_readonly("registeredNames", &Spheral::NodeListRegistrar<Dimension>::registeredNames)

    // Methods
    .def("valid", &Spheral::NodeListRegistrar<Dimension>::valid)
    ;

  //............................................................................
  // NodeList
  py::class_<NodeList<Dimension>,
             PyNodeList<NodeList<Dimension>>>(m, ("NodeList" + suffix).c_str())

    // Constructors
    .def(py::init<std::string, unsigned, unsigned, double, double, double, double, unsigned>(),
         "name"_a,
         "numInternal"_a = 0,
         "numGhost"_a = 0,
         "hmin"_a = 1.0e-20,
         "hmax"_a = 1.0e20,
         "hminratio"_a = 0.1,
         "nPerh"_a = 2.01,
         "maxNumNeighbors"_a = 500)

    // Methods
    .def("nodeType", &NodeList<Dimension>::nodeType)

    .def("mass", (Field<Dimension, Scalar>& (NodeList<Dimension>::*)()) &NodeList<Dimension>::mass, py::return_value_policy::reference_internal)
    .def("mass", (void (NodeList<Dimension>::*)(const Field<Dimension, Scalar>&)) &NodeList<Dimension>::mass)

    .def("positions", (Field<Dimension, Vector>& (NodeList<Dimension>::*)()) &NodeList<Dimension>::positions, py::return_value_policy::reference_internal)
    .def("positions", (void (NodeList<Dimension>::*)(const Field<Dimension, Vector>&)) &NodeList<Dimension>::positions)

    .def("velocity", (Field<Dimension, Vector>& (NodeList<Dimension>::*)()) &NodeList<Dimension>::velocity, py::return_value_policy::reference_internal)
    .def("velocity", (void (NodeList<Dimension>::*)(const Field<Dimension, Vector>&)) &NodeList<Dimension>::velocity)

    .def("Hfield", (Field<Dimension, SymTensor>& (NodeList<Dimension>::*)()) &NodeList<Dimension>::Hfield, py::return_value_policy::reference_internal)
    .def("Hfield", (void (NodeList<Dimension>::*)(const Field<Dimension, SymTensor>&)) &NodeList<Dimension>::Hfield)

    .def("work", (Field<Dimension, Scalar>& (NodeList<Dimension>::*)() const) &NodeList<Dimension>::work, py::return_value_policy::reference_internal)
    .def("work", (void (NodeList<Dimension>::*)(const Field<Dimension, Scalar>&)) &NodeList<Dimension>::work)

    .def("Hinverse", &NodeList<Dimension>::Hinverse)

    .def("neighbor", &NodeList<Dimension>::neighbor, py::return_value_policy::reference_internal)
    .def("registerNeighbor", &NodeList<Dimension>::registerNeighbor)
    .def("unregisterNeighbor", &NodeList<Dimension>::unregisterNeighbor)

    .def("haveField", &NodeList<Dimension>::haveField)

    // Virtual methods
    .def("deleteNodes", &NodeList<Dimension>::deleteNodes)
    .def("reorderNodes", &NodeList<Dimension>::reorderNodes)
    .def("label", &NodeList<Dimension>::label)
    .def("dumpState", &NodeList<Dimension>::dumpState)
    .def("restoreState", &NodeList<Dimension>::restoreState)

    // Attributes
    .def_property_readonly("name", &NodeList<Dimension>::name)
    .def_property_readonly("numNodes", &NodeList<Dimension>::numNodes)
    .def_property("numInternalNodes",
                  (unsigned (NodeList<Dimension>::*)() const) &NodeList<Dimension>::numInternalNodes,
                  (void (NodeList<Dimension>::*)(unsigned)) &NodeList<Dimension>::numInternalNodes)
    .def_property("numGhostNodes",
                  (unsigned (NodeList<Dimension>::*)() const) &NodeList<Dimension>::numGhostNodes,
                  (void (NodeList<Dimension>::*)(unsigned)) &NodeList<Dimension>::numGhostNodes)
    .def_property_readonly("numFields", &NodeList<Dimension>::numFields)
    .def_property_readonly("firstGhostNode", &NodeList<Dimension>::firstGhostNode)
    .def_property("nodesPerSmoothingScale",
                  (double (NodeList<Dimension>::*)() const) &NodeList<Dimension>::nodesPerSmoothingScale,
                  (void (NodeList<Dimension>::*)(double)) &NodeList<Dimension>::nodesPerSmoothingScale)
    .def_property("maxNumNeighbors",
                  (unsigned (NodeList<Dimension>::*)() const) &NodeList<Dimension>::maxNumNeighbors,
                  (void (NodeList<Dimension>::*)(unsigned)) &NodeList<Dimension>::maxNumNeighbors)
    .def_property("hmin",
                  (double (NodeList<Dimension>::*)() const) &NodeList<Dimension>::hmin,
                  (void (NodeList<Dimension>::*)(double)) &NodeList<Dimension>::hmin)
    .def_property("hmax",
                  (double (NodeList<Dimension>::*)() const) &NodeList<Dimension>::hmax,
                  (void (NodeList<Dimension>::*)(double)) &NodeList<Dimension>::hmax)
    .def_property("hminratio",
                  (double (NodeList<Dimension>::*)() const) &NodeList<Dimension>::hminratio,
                  (void (NodeList<Dimension>::*)(double)) &NodeList<Dimension>::hminratio)

    // Comparisons
    .def(py::self == py::self)
    .def(py::self != py::self)
    ;

  //............................................................................
  // FluidNodeList
  py::class_<FluidNodeList<Dimension>,
             PyFluidNodeList<Dimension, FluidNodeList<Dimension>>>(m, ("FluidNodeList" + suffix).c_str())

    // Constructors
    .def(py::init<std::string, EquationOfState<Dimension>&, int, int, double, double, double, double, int, double, double>(),
         "name"_a,
         "eos"_a,
         "numInternal"_a = 0,
         "numGhost"_a = 0,
         "hmin"_a = 1.0e-20,
         "hmax"_a = 1.0e20,
         "hminratio"_a = 0.1,
         "nPerh"_a = 2.01,
         "maxNumNeighbors"_a = 500,
         "rhoMin"_a  = 1.0e-10,
         "rhoMax"_a = 1.0e100)

    // Methods.
    .def("massDensity", (Field<Dimension, Scalar>& (FluidNodeList<Dimension>::*)()) &FluidNodeList<Dimension>::massDensity, py::return_value_policy::reference_internal)
    .def("massDensity", (void (FluidNodeList<Dimension>::*)(const Field<Dimension, Scalar>&)) &FluidNodeList<Dimension>::massDensity)
    
    .def("specificThermalEnergy", (Field<Dimension, Scalar>& (FluidNodeList<Dimension>::*)()) &FluidNodeList<Dimension>::specificThermalEnergy, py::return_value_policy::reference_internal)
    .def("specificThermalEnergy", (void (FluidNodeList<Dimension>::*)(const Field<Dimension, Scalar>&)) &FluidNodeList<Dimension>::specificThermalEnergy)
    
    .def("equationOfState", (const EquationOfState<Dimension>& (FluidNodeList<Dimension>::*)() const) &FluidNodeList<Dimension>::equationOfState)
    .def("equationOfState", (void (FluidNodeList<Dimension>::*)(const EquationOfState<Dimension>&)) &FluidNodeList<Dimension>::equationOfState)

    // Virtual methods
    .def("pressure", &FluidNodeList<Dimension>::pressure)
    .def("temperature", &FluidNodeList<Dimension>::temperature)
    .def("soundSpeed", &FluidNodeList<Dimension>::soundSpeed)
    .def("volume", &FluidNodeList<Dimension>::volume)
    .def("linearMomentum", &FluidNodeList<Dimension>::linearMomentum)
    .def("totalEnergy", &FluidNodeList<Dimension>::totalEnergy)

    // Attributes
    .def_property("rhoMin",
                  (double (FluidNodeList<Dimension>::*)() const) &FluidNodeList<Dimension>::rhoMin,
                  (double (FluidNodeList<Dimension>::*)() const) &FluidNodeList<Dimension>::rhoMin)
    .def_property("rhoMax",
                  (double (FluidNodeList<Dimension>::*)() const) &FluidNodeList<Dimension>::rhoMax,
                  (double (FluidNodeList<Dimension>::*)() const) &FluidNodeList<Dimension>::rhoMax)
    ;

  //............................................................................
  // SmoothingScaleBase
  py::class_<SmoothingScaleBase<Dimension>,
             PySmoothingScale<Dimension, SmoothingScaleBase<Dimension>>> SmoothingScaleBasePB11(m, ("SmoothingScaleBase" + suffix).c_str());
  smoothingScaleObjectBindings<Dimension, SmoothingScaleBase<Dimension>>(m, suffix, SmoothingScaleBasePB11);

  //............................................................................
  // FixedSmoothingScale
  py::class_<FixedSmoothingScale<Dimension>, SmoothingScaleBase<Dimension>, 
             PySmoothingScale<Dimension, FixedSmoothingScale<Dimension>>> FixedSmoothingScalePB11(m, ("FixedSmoothingScale" + suffix).c_str());
  smoothingScaleObjectBindings<Dimension, FixedSmoothingScale<Dimension>>(m, suffix, FixedSmoothingScalePB11);

  //............................................................................
  // SPHSmoothingScale
  py::class_<SPHSmoothingScale<Dimension>, SmoothingScaleBase<Dimension>, 
             PySmoothingScale<Dimension, SPHSmoothingScale<Dimension>>> SPHSmoothingScalePB11(m, ("SPHSmoothingScale" + suffix).c_str());
  smoothingScaleObjectBindings<Dimension, SPHSmoothingScale<Dimension>>(m, suffix, SPHSmoothingScalePB11);

  //............................................................................
  // ASPHSmoothingScale
  py::class_<ASPHSmoothingScale<Dimension>, SmoothingScaleBase<Dimension>, 
             PySmoothingScale<Dimension, ASPHSmoothingScale<Dimension>>> ASPHSmoothingScalePB11(m, ("ASPHSmoothingScale" + suffix).c_str());
  smoothingScaleObjectBindings<Dimension, ASPHSmoothingScale<Dimension>>(m, suffix, ASPHSmoothingScalePB11);

  //............................................................................
  // STL containers
  py::bind_vector<std::vector<NodeList<Dimension>*>>(m, "vector_of_NodeList" + suffix);
  py::bind_vector<std::vector<FluidNodeList<Dimension>*>>(m, "vector_of_FluidNodeList" + suffix);

  //............................................................................
  // Methods
  m.def("generateVoidNodes", &generateVoidNodes<Dimension>);
  m.def("zerothNodalMoment", &nthNodalMoment<Dimension, 0U>);
  m.def("firstNodalMoment", &nthNodalMoment<Dimension, 1>);
  m.def("secondNodalMoment", &nthNodalMoment<Dimension, 2U>);
  m.def("zerothAndFirstNodalMoments", &zerothAndFirstNodalMoments<Dimension>);

}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralNodeList) {
  py::module m("SpheralNodeList", "Spheral NodeList module.");

  //............................................................................
  // NodeType
  py::enum_<Spheral::NodeSpace::NodeType>(m, "NodeType")
    .value("InternalNode", Spheral::NodeSpace::NodeType::InternalNode)
    .value("GhostNode", Spheral::NodeSpace::NodeType::GhostNode)
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
