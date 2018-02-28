// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "Physics/GenericBodyForce.hh"
#include "ExternalForce/PointPotential.hh"
#include "ExternalForce/ConstantAcceleration.hh"
#include "ExternalForce/LinearAcceleration.hh"

#include "Pybind11Wraps/Physics/PyGenericBodyForce.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::PhysicsSpace;

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

  //............................................................................
  // PointPotential
  typedef PointPotential<Dimension> PP;
  py::class_<PP,
             PyGenericBodyForce<Dimension, GenericBodyForce<Dimension> > > ppPB11(m, ("PointPotential" + suffix).c_str());
    
  // Constructors
  ppPB11.def(py::init<double, double, double, const Vector&>(), "G"_a, "mass"_a, "coreRadius"_a, "origin"_a);

  // Methods
  ppPB11.def("specificPotential", &PP::specificPotential, "r"_a);

  // Attributes
  ppPB11.def_property("G", &PP::G, &PP::setG);
  ppPB11.def_property("mass", &PP::mass, &PP::setMass);
  ppPB11.def_property("coreRadius", &PP::coreRadius, &PP::setCoreRadius);
  ppPB11.def_property("origin", &PP::origin, &PP::setOrigin);
  ppPB11.def_property("deltaPotentialFraction", &PP::deltaPotentialFraction, &PP::setDeltaPotentialFraction);
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralExternalForce, m) {
  using namespace Spheral;
  using namespace Spheral::PhysicsSpace;

  m.doc() = "Spheral ExternalForce module.";

  //............................................................................
  // Per dimension bindings.
#ifdef SPHERAL1D
  dimensionBindings<Spheral::Dim<1>>(m, "1d");
#endif

// #ifdef SPHERAL2D
//   dimensionBindings<Spheral::Dim<2>>(m, "2d");
// #endif

// #ifdef SPHERAL3D
//   dimensionBindings<Spheral::Dim<3>>(m, "3d");
// #endif
}
