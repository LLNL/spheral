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
#include "PyGenericHydro.hh"
#include "PyGenericBodyForce.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"
#include "virtualPhysicsBindings.hh"

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
    
  // Constructors
  phyPB11.def(py::init<>());

  // Methods
  phyPB11.def("appendBoundary", &Phys::appendBoundary, "boundary"_a);
  phyPB11.def("prependBoundary", &Phys::prependBoundary, "boundary"_a);
  phyPB11.def("clearBoundaries", &Phys::clearBoundaries);
  phyPB11.def("haveBoundary", &Phys::haveBoundary, "boundary"_a);
  phyPB11.def("boundaryConditions", &Phys::boundaryConditions);

  //............................................................................
  // GenericHydro
  typedef GenericHydro<Dimension> GH;
  py::class_<GH, Phys,
             PyGenericHydro<Dimension, GH>> ghPB11(m, ("GenericHydro" + suffix).c_str());

  // Constructors
  ghPB11.def(py::init<const TableKernel<Dimension>&, const TableKernel<Dimension>&, ArtificialViscosity<Dimension>&, const double, const bool>(),
             "W"_a, "WPi"_a, "Q"_a, "cfl"_a, "useVelocityMagnitudeForDt"_a);

  // Methods
  ghPB11.def("artificialViscosity", &GH::artificialViscosity);
  ghPB11.def("kernel", &GH::kernel);
  ghPB11.def("PiKernel", &GH::PiKernel);

  // Virtual methods
  ghPB11.def("dt", &GH::dt, "dataBase"_a, "state"_a, "derivs"_a, "currentTime"_a);

  // Attributes
  ghPB11.def_property("cfl", (Scalar (GH::*)() const) &GH::cfl, (void (GH::*)(Scalar)) &GH::cfl);
  ghPB11.def_property("useVelocityMagnitudeForDt", (bool (GH::*)() const) &GH::useVelocityMagnitudeForDt, (void (GH::*)(bool)) &GH::useVelocityMagnitudeForDt);
  ghPB11.def_property_readonly("minMasterNeighbor", &GH::minMasterNeighbor);
  ghPB11.def_property_readonly("maxMasterNeighbor", &GH::maxMasterNeighbor);
  ghPB11.def_property_readonly("averageMasterNeighbor", &GH::averageMasterNeighbor);
  ghPB11.def_property_readonly("minCoarseNeighbor", &GH::minCoarseNeighbor);
  ghPB11.def_property_readonly("maxCoarseNeighbor", &GH::maxCoarseNeighbor);
  ghPB11.def_property_readonly("averageCoarseNeighbor", &GH::averageCoarseNeighbor);
  ghPB11.def_property_readonly("minRefineNeighbor", &GH::minRefineNeighbor);
  ghPB11.def_property_readonly("maxRefineNeighbor", &GH::maxRefineNeighbor);
  ghPB11.def_property_readonly("averageRefineNeighbor", &GH::averageRefineNeighbor);
  ghPB11.def_property_readonly("minActualNeighbor", &GH::minActualNeighbor);
  ghPB11.def_property_readonly("maxActualNeighbor", &GH::maxActualNeighbor);
  ghPB11.def_property_readonly("averageActualNeighbor", &GH::averageActualNeighbor);

  //............................................................................
  // GenericBodyForce
  typedef GenericBodyForce<Dimension> GBF;
  py::class_<GBF, Phys,
             PyGenericBodyForce<Dimension, GBF> > gbfPB11(m, ("GenericBodyForce" + suffix).c_str());

  // Constructors
  gbfPB11.def(py::init<>());

  // Methods
  gbfPB11.def("DxDt", &GBF::DxDt);
  gbfPB11.def("DvDt", &GBF::DvDt);

  // Virtual methods
  gbfPB11.def("registerState", &GBF::registerState);
  gbfPB11.def("registerDerivatives", &GBF::registerDerivatives);

  //............................................................................
  // The STL containers of Physics objects.
  py::bind_vector<std::vector<Phys*>>(m, "vector_of_Physics" + suffix);
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralPhysics, m) {
  using namespace Spheral;
  using namespace Spheral::PhysicsSpace;

  m.doc() = "Spheral Physics module.";

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
