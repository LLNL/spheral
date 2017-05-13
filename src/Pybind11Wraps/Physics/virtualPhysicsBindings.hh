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

  obj

    // Methods
    .def("evaluateDerivatives", &Obj::evaluateDerivatives, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a)
    .def("dt", &Obj::dt, "dataBase"_a, "state"_a, "derivs"_a, "currentTime"_a)
    .def("registerState", &Obj::registerState, "dataBase"_a, "state"_a)
    .def("registerDerivatives", &Obj::registerDerivatives, "dataBase"_a, "state"_a)
    .def("label", &Obj::label)
    .def("applyGhostBoundaries", &Obj::applyGhostBoundaries, "state"_a, "derivs"_a)
    .def("enforceBoundaries", &Obj::enforceBoundaries, "state"_a, "derivs"_a)
    .def("initializeProblemStartup", &Obj::initializeProblemStartup)
    .def("preStepInitialize", &Obj::preStepInitialize, "dataBase"_a, "state"_a, "derivs"_a)
    .def("initialize", &Obj::initialize, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a)
    .def("finalize", &Obj::finalize, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a)
    .def("finalizeDerivatives", &Obj::finalizeDerivatives, "time"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a)
    .def("postStateUpdate", &Obj::postStateUpdate, "dataBase"_a, "state"_a, "derivs"_a)
    .def("requireConnectivity", &Obj::requireConnectivity)
    .def("requireGhostConnectivity", &Obj::requireGhostConnectivity)
    .def("extraEnergy", &Obj::extraEnergy)
    .def("extraMomentum", &Obj::extraMomentum)
    ;
}

}
}

#endif
