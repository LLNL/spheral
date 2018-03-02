// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "Physics/GenericBodyForce.hh"
#include "Gravity/NBodyGravity.hh"
#include "Gravity/TreeGravity.hh"

#include "Pybind11Wraps/Physics/PyGenericBodyForce.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::PhysicsSpace;
using namespace Spheral::GravitySpace;

namespace {  // anonymous

//------------------------------------------------------------------------------
// NBodyGravity
//------------------------------------------------------------------------------
template<typename Dimension>
void NBodyGravityBindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  typedef NBodyGravity<Dimension> NBG;
  py::class_<NBG, PyGenericBodyForce<Dimension, NBG> > nbgPB11(m, ("NBodyGravity" + suffix).c_str());
    
  // Constructors
  nbgPB11.def(py::init<const double, const double, const double, const bool>(),
              "plummerSofteningLength"_a, "maxDeltaVelocity"_a, "G"_a, "compatibleVelocityUpdate"_a);

  // Attributes
  nbgPB11.def_property_readonly("G", &NBG::G);
  nbgPB11.def_property("softeningLength",
                       (double (NBG::*)() const) &NBG::softeningLength,
                       (void (NBG::*)(const double)) &NBG::softeningLength);
  nbgPB11.def_property("compatibleVelocityUpdate",
                       (bool (NBG::*)() const) &NBG::compatibleVelocityUpdate,
                       (void (NBG::*)(const bool)) &NBG::compatibleVelocityUpdate);
}

//------------------------------------------------------------------------------
// TreeGravity
//------------------------------------------------------------------------------
template<typename Dimension>
void TreeGravityBindings(py::module& m, const std::string prefix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  typedef TreeGravity<Dimension> TG;
  py::class_<TG, PyGenericBodyForce<Dimension, TG> > tgPB11(m, (prefix + "TreeGravity").c_str());
    
  // Constructors
  tgPB11.def(py::init<const double, const double, const double, const double, const GravityTimeStepType>(),
             "G"_a, "softeningLength"_a, "opening"_a, "ftimestep"_a, "timeStepChoice"_a);

  // Attributes
  tgPB11.def_property_readonly("G", &TG::G);
  tgPB11.def_property("opening",
                      (double (TG::*)() const) &TG::opening,
                      (void (TG::*)(const double)) &TG::opening);
  tgPB11.def_property("softeningLength",
                      (double (TG::*)() const) &TG::softeningLength,
                      (void (TG::*)(const double)) &TG::softeningLength);
  tgPB11.def_property("fimestep",
                      (double (TG::*)() const) &TG::ftimestep,
                      (void (TG::*)(const double)) &TG::ftimestep);
  tgPB11.def_property("timeStepChoice",
                      (GravityTimeStepType (TG::*)() const) &TG::timeStepChoice,
                      (void (TG::*)(const GravityTimeStepType)) &TG::timeStepChoice);
  tgPB11.def_property_readonly("xmin", &TG::xmin);
  tgPB11.def_property_readonly("xmax", &TG::xmax);
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralGravity, m) {
  using namespace Spheral;
  using namespace Spheral::PhysicsSpace;

  m.doc() = "Spheral Gravity module.";

  //............................................................................
  // GravityType
  py::enum_<Spheral::GravitySpace::GravityTimeStepType>(m, "GravityTimeStepType")
    .value("AccelerationRatio", Spheral::GravitySpace::GravityTimeStepType::AccelerationRatio)
    .value("DynamicalTime", Spheral::GravitySpace::GravityTimeStepType::DynamicalTime)
    .export_values();

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  NBodyGravityBindings<Spheral::Dim<1>>(m, "1d");
#endif

#ifdef SPHERAL2D
  NBodyGravityBindings<Spheral::Dim<2>>(m, "2d");
  TreeGravityBindings<Spheral::Dim<2>>(m, "Quad");
#endif

#ifdef SPHERAL3D
  NBodyGravityBindings<Spheral::Dim<3>>(m, "3d");
  TreeGravityBindings<Spheral::Dim<3>>(m, "Oct");
#endif
}
