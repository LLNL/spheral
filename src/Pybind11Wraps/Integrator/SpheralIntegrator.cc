// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Physics/Physics.hh"
#include "Integrator/Integrator.hh"
#include "Integrator/PredictorCorrector.hh"
#include "Integrator/SynchronousRK1.hh"
#include "Integrator/SynchronousRK2.hh"
#include "Integrator/SynchronousRK4.hh"
#include "Integrator/CheapSynchronousRK2.hh"
#include "Integrator/Verlet.hh"

#include "Pybind11Wraps/Integrator/PyAbstractIntegrator.hh"
#include "Pybind11Wraps/Integrator/PyIntegrator.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::PhysicsSpace;
using namespace Spheral::IntegratorSpace;
using DataBaseSpace::DataBase;
using PhysicsSpace::Physics;

namespace {  // anonymous

//------------------------------------------------------------------------------
// Generic integrator methods.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void integratorMethodsBindings(PB11Obj& obj) {
  typedef typename Dimension::Scalar Scalar;
    
  // Constructors
  obj.def(py::init<>());
  obj.def(py::init<DataBase<Dimension>&>(), "dataBase"_a);
  obj.def(py::init<DataBase<Dimension>&, const std::vector<Physics<Dimension>*>&>(),
          "dataBase"_a, "physicsPackages"_a);

  // Methods
  obj.def("step", (void (Obj::*)(Scalar, State<Dimension>&, StateDerivatives<Dimension>&)) &Obj::step,
          "maxTime"_a, "state"_a, "derivs"_a);
}

//------------------------------------------------------------------------------
// Bind the objects templated on Dimension.
//------------------------------------------------------------------------------
template<typename Dimension>
void dimensionBindings(py::module& m, const std::string suffix) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  //............................................................................
  typedef Integrator<Dimension> INT;
  py::class_<INT, PyRestartMethods<PyAbstractIntegrator<Dimension, INT>>>  intPB11(m, ("Integrator" + suffix).c_str());
  restartMethodBindings<INT>(intPB11);   // Bind restart methods.
    
  // Constructors
  intPB11.def(py::init<>());
  intPB11.def(py::init<DataBase<Dimension>&>(), "dataBase"_a);
  intPB11.def(py::init<DataBase<Dimension>&, const std::vector<Physics<Dimension>*>&>(),
              "dataBase"_a, "physicsPackages"_a);

  // Methods
  intPB11.def("step", (void (INT::*)(Scalar, State<Dimension>&, StateDerivatives<Dimension>&)) &INT::step,
              "maxTime"_a, "state"_a, "derivs"_a);
  intPB11.def("step", (void (INT::*)(Scalar)) &INT::step, "maxTime"_a);
  intPB11.def("selectDt", &INT::selectDt, "dtMin"_a, "dtMax"_a, "state"_a, "derivs"_a);
  intPB11.def("preStepInitialize", &INT::preStepInitialize, "state"_a, "derivs"_a);
  intPB11.def("initializeDerivatives", &INT::initializeDerivatives, "t"_a, "dt"_a, "state"_a, "derivs"_a);
  intPB11.def("evaluateDerivatives", &INT::evaluateDerivatives, "t"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
  intPB11.def("finalizeDerivatives", &INT::finalizeDerivatives, "t"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
  intPB11.def("postStateUpdate", &INT::postStateUpdate, "dataBase"_a, "state"_a, "derivs"_a);
  intPB11.def("finalizeDerivatives", &INT::finalizeDerivatives, "t"_a, "dt"_a, "dateBase"_a, "state"_a, "derivs"_a);
  intPB11.def("appendPhysicsPackage", &INT::appendPhysicsPackage, "package"_a);
  intPB11.def("havePhysicsPackage", &INT::havePhysicsPackage, "package"_a);
  intPB11.def("uniqueBoundaryConditions", &INT::uniqueBoundaryConditions);
  intPB11.def("setGhostNodes", &INT::setGhostNodes);
  intPB11.def("applyGhostBoundaries", &INT::applyGhostBoundaries, "state"_a, "derivs"_a);
  intPB11.def("finalizeGhostBoundaries", &INT::finalizeGhostBoundaries);
  intPB11.def("setViolationNodes", &INT::setViolationNodes);
  intPB11.def("enforceBoundaries", &INT::enforceBoundaries, "state"_a, "derivs"_a);
  intPB11.def("copyGhostState", &INT::copyGhostState, "state0"_a, "state1"_a);
  intPB11.def("advance", &INT::advance, "goalTime"_a);

  // Attributes
  intPB11.def_property("currentTime",
                       (Scalar (INT::*)() const) &INT::currentTime,
                       (void (INT::*)(Scalar)) &INT::currentTime);
  intPB11.def_property("currentCycle",
                       (int (INT::*)() const) &INT::currentCycle,
                       (void (INT::*)(int)) &INT::currentCycle);
  intPB11.def_property("dtMin",
                       (Scalar (INT::*)() const) &INT::dtMin,
                       (void (INT::*)(Scalar)) &INT::dtMin);
  intPB11.def_property("dtMax",
                       (Scalar (INT::*)() const) &INT::dtMax,
                       (void (INT::*)(Scalar)) &INT::dtMax);
  intPB11.def_property("lastDt",
                       (Scalar (INT::*)() const) &INT::lastDt,
                       (void (INT::*)(Scalar)) &INT::lastDt);
  intPB11.def_property("dtGrowth",
                       (Scalar (INT::*)() const) &INT::dtGrowth,
                       (void (INT::*)(Scalar)) &INT::dtGrowth);
  intPB11.def_property_readonly("dataBase", &INT::dataBase);
  intPB11.def_property_readonly("physicsPackages", &INT::physicsPackages);
  intPB11.def_property("rigorousBoundaries",
                       (bool (INT::*)() const) &INT::rigorousBoundaries,
                       (void (INT::*)(bool)) &INT::rigorousBoundaries);
  intPB11.def_property("updateBoundaryFrequency",
                       (int (INT::*)() const) &INT::updateBoundaryFrequency,
                       (void (INT::*)(const int)) &INT::updateBoundaryFrequency);
  intPB11.def_property("verbose",
                       (bool (INT::*)() const) &INT::verbose,
                       (void (INT::*)(bool)) &INT::verbose);
  intPB11.def_property("domainDecompositionIndependent",
                       (bool (INT::*)() const) &INT::domainDecompositionIndependent,
                       (void (INT::*)(bool)) &INT::domainDecompositionIndependent);
  intPB11.def_property("cullGhostNodes",
                       (bool (INT::*)() const) &INT::cullGhostNodes,
                       (void (INT::*)(bool)) &INT::cullGhostNodes);

  // // Protected methods
  // intPB11.def("accessDataBase", &INT::accessDataBase);

  //............................................................................
  typedef PredictorCorrector<Dimension> PC;
  py::class_<PC, PyRestartMethods<PyIntegrator<Dimension, PC>>>  pcPB11(m, ("PredictorCorrectorIntegrator" + suffix).c_str());
  integratorMethodsBindings<Dimension, PC>(pcPB11);
    
  //............................................................................
  typedef SynchronousRK1<Dimension> SRK1;
  py::class_<SRK1, PyRestartMethods<PyIntegrator<Dimension, SRK1>>>  srk1PB11(m, ("SynchronousRK1Integrator" + suffix).c_str());
  integratorMethodsBindings<Dimension, SRK1>(srk1PB11);
    
  //............................................................................
  typedef SynchronousRK2<Dimension> SRK2;
  py::class_<SRK2, PyRestartMethods<PyIntegrator<Dimension, SRK2>>>  srk2PB11(m, ("SynchronousRK2Integrator" + suffix).c_str());
  integratorMethodsBindings<Dimension, SRK2>(srk2PB11);
    
  //............................................................................
  typedef SynchronousRK4<Dimension> SRK4;
  py::class_<SRK4, PyRestartMethods<PyIntegrator<Dimension, SRK4>>>  srk4PB11(m, ("SynchronousRK4Integrator" + suffix).c_str());
  integratorMethodsBindings<Dimension, SRK4>(srk4PB11);
    
  //............................................................................
  typedef CheapSynchronousRK2<Dimension> CSRK2;
  py::class_<CSRK2, PyRestartMethods<PyIntegrator<Dimension, CSRK2>>>  csrk2PB11(m, ("CheapSynchronousRK2Integrator" + suffix).c_str());
  integratorMethodsBindings<Dimension, CSRK2>(csrk2PB11);
    
  //............................................................................
  typedef Verlet<Dimension> VRL;
  py::class_<VRL, PyRestartMethods<PyIntegrator<Dimension, VRL>>>  vrlPB11(m, ("VerletIntegrator" + suffix).c_str());
  integratorMethodsBindings<Dimension, VRL>(vrlPB11);
    
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralGravity, m) {
  using namespace Spheral;
  using namespace Spheral::IntegratorSpace;

  m.doc() = "Spheral Integrator module: provides time integration methods.";

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
