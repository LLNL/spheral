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

#include "PyAbstractPhysics.hh"
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

  //............................................................................
  // Physics
  typedef Physics<Dimension> Phy;
  py::class_<Phy,
             PyAbstractPhysics<Dimension, Phy>> phyPB11(m, ("Physics" + suffix).c_str());
  virtualPhysicsBindings<Dimension, Phy>(m, phyPB11);
  phyPB11
    
    // Constructors
    .def(py::init<>())
    ;

  //............................................................................
  // The STL containers of Physics objects.
  py::bind_vector<std::vector<Phy*>>(m, "vector_of_Physics" + suffix);
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
