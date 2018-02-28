#ifndef __Spheral_virtualPhysicsBindings__
#define __Spheral_virtualPhysicsBindings__

namespace Spheral {
namespace PhysicsSpace {

//------------------------------------------------------------------------------
// Common virtual methods of Physics objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void virtualPhysicsBindings(py::module& m, PB11Obj& obj) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Methods
  obj.def("evaluateDerivatives", &Obj::evaluateDerivatives, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
  obj.def("dt", &Obj::dt, "dataBase"_a, "state"_a, "derivs"_a, "currentTime"_a);
  obj.def("registerState", &Obj::registerState, "dataBase"_a, "state"_a);
  obj.def("registerDerivatives", &Obj::registerDerivatives, "dataBase"_a, "state"_a);
  obj.def("label", &Obj::label);
  obj.def("applyGhostBoundaries", &Obj::applyGhostBoundaries, "state"_a, "derivs"_a);
  obj.def("enforceBoundaries", &Obj::enforceBoundaries, "state"_a, "derivs"_a);
  obj.def("initializeProblemStartup", &Obj::initializeProblemStartup);
  obj.def("preStepInitialize", &Obj::preStepInitialize, "dataBase"_a, "state"_a, "derivs"_a);
  obj.def("initialize", &Obj::initialize, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
  obj.def("finalize", &Obj::finalize, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
  obj.def("finalizeDerivatives", &Obj::finalizeDerivatives, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
  obj.def("postStateUpdate", &Obj::postStateUpdate, "dataBase"_a, "state"_a, "derivs"_a);
  obj.def("requireConnectivity", &Obj::requireConnectivity);
  obj.def("requireGhostConnectivity", &Obj::requireGhostConnectivity);
  obj.def("extraEnergy", &Obj::extraEnergy);
  obj.def("extraMomentum", &Obj::extraMomentum);
}

}
}

#endif
