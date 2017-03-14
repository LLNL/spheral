// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "Physics/Physics.hh"
#include "Physics/GenericHydro.hh"
#include "Physics/GenericBodyForce.hh"
#include "Boundary/Boundary.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "Kernel/TableKernel.hh"

#include "PyAbstractPhysics.hh"
#include "PyPhysics.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::PhysicsSpace;

//------------------------------------------------------------------------------
// 1D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Physics<Dim<1>>*>);

//------------------------------------------------------------------------------
// 2D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Physics<Dim<2>>*>);

//------------------------------------------------------------------------------
// 3D
//------------------------------------------------------------------------------
PYBIND11_MAKE_OPAQUE(std::vector<Physics<Dim<3>>*>);

namespace {  // anonymous

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
  typedef Spheral::GeomPlane<Dimension> Plane;
  using Spheral::NodeSpace::NodeList;
  using Spheral::KernelSpace::TableKernel;
  using Spheral::ArtificialViscositySpace::ArtificialViscosity;

  //............................................................................
  // Physics
  typedef Physics<Dimension> Phys;
  py::class_<Phys,
             PyAbstractPhysics<Dimension, Phys>> phyPB11(m, ("Physics" + suffix).c_str());
  virtualPhysicsBindings<Dimension, Phys>(m, phyPB11);
  phyPB11
    
    // Constructors
    .def(py::init<>())

    // Methods
    .def("appendBoundary", &Phys::appendBoundary, "boundary"_a)
    .def("prependBoundary", &Phys::prependBoundary, "boundary"_a)
    .def("clearBoundaries", &Phys::clearBoundaries)
    .def("haveBoundary", &Phys::haveBoundary, "boundary"_a)
    .def("boundaryConditions", &Phys::boundaryConditions)
    ;

  //............................................................................
  // GenericHydro
  typedef GenericHydro<Dimension> GH;
  py::class_<GH, Phys,
             PyAbstractPhysics<Dimension, GH>> ghPB11(m, ("GenericHydro" + suffix).c_str());
  //  virtualPhysicsBindings<Dimension, Phys>(m, ghPB11);
  ghPB11

    // Constructors
    .def(py::init<const TableKernel<Dimension>&, const TableKernel<Dimension>&, ArtificialViscosity<Dimension>&, const double, const bool>(),
         "W"_a, "WPi"_a, "Q"_a, "cfl"_a, "useVelocityMagnitudeForDt"_a)

    // Methods
    .def("artificialViscosity", &GH::artificialViscosity)
    .def("kernel", &GH::kernel)
    .def("PiKernel", &GH::PiKernel)

    // Virtual methods
    .def("dt", &GH::dt, "dataBase"_a, "state"_a, "derivs"_a, "currentTime"_a)

    // Attributes
    .def_property("cfl", (Scalar (GH::*)() const) &GH::cfl, (void (GH::*)(Scalar)) &GH::cfl)
    .def_property("useVelocityMagnitudeForDt", (bool (GH::*)() const) &GH::useVelocityMagnitudeForDt, (void (GH::*)(bool)) &GH::useVelocityMagnitudeForDt)
    .def_property_readonly("minMasterNeighbor", &GH::minMasterNeighbor)
    .def_property_readonly("maxMasterNeighbor", &GH::maxMasterNeighbor)
    .def_property_readonly("averageMasterNeighbor", &GH::averageMasterNeighbor)
    .def_property_readonly("minCoarseNeighbor", &GH::minCoarseNeighbor)
    .def_property_readonly("maxCoarseNeighbor", &GH::maxCoarseNeighbor)
    .def_property_readonly("averageCoarseNeighbor", &GH::averageCoarseNeighbor)
    .def_property_readonly("minRefineNeighbor", &GH::minRefineNeighbor)
    .def_property_readonly("maxRefineNeighbor", &GH::maxRefineNeighbor)
    .def_property_readonly("averageRefineNeighbor", &GH::averageRefineNeighbor)
    .def_property_readonly("minActualNeighbor", &GH::minActualNeighbor)
    .def_property_readonly("maxActualNeighbor", &GH::maxActualNeighbor)
    .def_property_readonly("averageActualNeighbor", &GH::averageActualNeighbor)
    ;

  //............................................................................
  // The STL containers of Physics objects.
  py::bind_vector<std::vector<Phys*>>(m, "vector_of_Physics" + suffix);
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralPhysics) {
  using namespace Spheral;
  using namespace Spheral::PhysicsSpace;

  py::module m("SpheralPhysics", "Spheral Physics module.");

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif

// #ifdef SPHERAL2D
//   dimensionBindings<Spheral::Dim<2>>(m, "2d");
//   twoDimensionalBindings(m);
// #endif

// #ifdef SPHERAL3D
//   dimensionBindings<Spheral::Dim<3>>(m, "3d");
//   threeDimensionalBindings(m);
// #endif

  return m.ptr();
}
