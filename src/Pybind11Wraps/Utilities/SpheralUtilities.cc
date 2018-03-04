// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>
#include <utility>

#include "Geometry/Dimension.hh"
#include "boost/math/special_functions/legendre.hpp"
#include "Utilities/Functors.hh"
#include "Utilities/erff.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/simpsonsIntegration.hh"
#include "Utilities/globalNodeIDs.hh"
#include "Utilities/rotationMatrix.hh"
#include "Utilities/iterateIdealH.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/peanoHilbertOrderIndices.hh"
#include "Utilities/boundingBox.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/testBoxIntersection.hh"
#include "Utilities/lineSegmentIntersections.hh"
#include "Utilities/pointDistances.hh"
#include "Utilities/segmentIntersectEdges.hh"
#include "Utilities/integrateThroughMeshAlongSegment.hh"
#include "Utilities/numberDensity.hh"
//#include "Utilities/writeRectilinearMesh.hh"
#include "Utilities/planarReflectingOperator.hh"
#include "Utilities/KeyTraits.hh"
#include "Utilities/pointInPolygon.hh"
#include "Utilities/pointOnPolygon.hh"
#include "Utilities/pointInPolyhedron.hh"
#include "Utilities/pointOnPolyhedron.hh"
#include "Utilities/refinePolyhedron.hh"
#include "Utilities/overlayRemapFields.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "Utilities/Timer.hh"

#include "Pybind11Wraps/Integrator/PyAbstractIntegrator.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

#ifndef NOR3D
#include "Utilities/r3d_utils.hh"
#else
//------------------------------------------------------------------------------
// Stub these methods out when we're not allowing R3D.
//------------------------------------------------------------------------------
namespace Spheral {
  inline Dim<2>::FacetedVolume clipFacetedVolume(const Dim<2>::FacetedVolume& poly,
                                                 const std::vector<GeomPlane<Dim<2> > >& planes) {
    VERIFY2(false, "ERROR: clipFacetedVolume unavailable without R3D.");
  }

  inline Dim<3>::FacetedVolume clipFacetedVolume(const Dim<3>::FacetedVolume& poly,
                                                 const std::vector<GeomPlane<Dim<3> > >& planes) {
    VERIFY2(false, "ERROR: clipFacetedVolume unavailable without R3D.");
  }
}
#endif

namespace py = pybind11;
using namespace pybind11::literals;
using Spheral::NodeSpace::NodeList;

using namespace Spheral;

namespace {  // anonymous

// //------------------------------------------------------------------------------
// // Generic integrator methods.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename Obj, typename PB11Obj>
// void integratorMethodsBindings(PB11Obj& obj) {
//   typedef typename Dimension::Scalar Scalar;
    
//   // Constructors
//   obj.def(py::init<>());
//   obj.def(py::init<DataBase<Dimension>&>(), "dataBase"_a);
//   obj.def(py::init<DataBase<Dimension>&, const std::vector<Physics<Dimension>*>&>(),
//           "dataBase"_a, "physicsPackages"_a);

//   // Methods
//   obj.def("step", (void (Obj::*)(Scalar, State<Dimension>&, StateDerivatives<Dimension>&)) &Obj::step,
//           "maxTime"_a, "state"_a, "derivs"_a);
// }

// //------------------------------------------------------------------------------
// // Bind the objects templated on Dimension.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void dimensionBindings(py::module& m, const std::string suffix) {

//   typedef typename Dimension::Scalar Scalar;
//   typedef typename Dimension::Vector Vector;
//   typedef typename Dimension::Tensor Tensor;
//   typedef typename Dimension::SymTensor SymTensor;
//   typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

//   //............................................................................
//   typedef Integrator<Dimension> INT;
//   py::class_<INT, PyRestartMethods<PyAbstractIntegrator<Dimension, INT>>>  intPB11(m, ("Integrator" + suffix).c_str());
//   restartMethodBindings<INT>(m, intPB11);   // Bind restart methods.
    
//   // Constructors
//   intPB11.def(py::init<>());
//   intPB11.def(py::init<DataBase<Dimension>&>(), "dataBase"_a);
//   intPB11.def(py::init<DataBase<Dimension>&, const std::vector<Physics<Dimension>*>&>(),
//               "dataBase"_a, "physicsPackages"_a);

//   // Methods
//   intPB11.def("step", (void (INT::*)(Scalar, State<Dimension>&, StateDerivatives<Dimension>&)) &INT::step,
//               "maxTime"_a, "state"_a, "derivs"_a);
//   intPB11.def("step", (void (INT::*)(Scalar)) &INT::step, "maxTime"_a);
//   intPB11.def("selectDt", &INT::selectDt, "dtMin"_a, "dtMax"_a, "state"_a, "derivs"_a);
//   intPB11.def("preStepInitialize", &INT::preStepInitialize, "state"_a, "derivs"_a);
//   intPB11.def("initializeDerivatives", &INT::initializeDerivatives, "t"_a, "dt"_a, "state"_a, "derivs"_a);
//   intPB11.def("evaluateDerivatives", &INT::evaluateDerivatives, "t"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
//   intPB11.def("finalizeDerivatives", &INT::finalizeDerivatives, "t"_a, "dt"_a, "dataBase"_a, "state"_a, "derivs"_a);
//   intPB11.def("postStateUpdate", &INT::postStateUpdate, "dataBase"_a, "state"_a, "derivs"_a);
//   intPB11.def("finalizeDerivatives", &INT::finalizeDerivatives, "t"_a, "dt"_a, "dateBase"_a, "state"_a, "derivs"_a);
//   intPB11.def("appendPhysicsPackage", &INT::appendPhysicsPackage, "package"_a);
//   intPB11.def("havePhysicsPackage", &INT::havePhysicsPackage, "package"_a);
//   intPB11.def("uniqueBoundaryConditions", &INT::uniqueBoundaryConditions);
//   intPB11.def("setGhostNodes", &INT::setGhostNodes);
//   intPB11.def("applyGhostBoundaries", &INT::applyGhostBoundaries, "state"_a, "derivs"_a);
//   intPB11.def("finalizeGhostBoundaries", &INT::finalizeGhostBoundaries);
//   intPB11.def("setViolationNodes", &INT::setViolationNodes);
//   intPB11.def("enforceBoundaries", &INT::enforceBoundaries, "state"_a, "derivs"_a);
//   intPB11.def("copyGhostState", &INT::copyGhostState, "state0"_a, "state1"_a);
//   intPB11.def("advance", &INT::advance, "goalTime"_a);

//   // Attributes
//   intPB11.def_property("currentTime",
//                        (Scalar (INT::*)() const) &INT::currentTime,
//                        (void (INT::*)(Scalar)) &INT::currentTime);
//   intPB11.def_property("currentCycle",
//                        (int (INT::*)() const) &INT::currentCycle,
//                        (void (INT::*)(int)) &INT::currentCycle);
//   intPB11.def_property("dtMin",
//                        (Scalar (INT::*)() const) &INT::dtMin,
//                        (void (INT::*)(Scalar)) &INT::dtMin);
//   intPB11.def_property("dtMax",
//                        (Scalar (INT::*)() const) &INT::dtMax,
//                        (void (INT::*)(Scalar)) &INT::dtMax);
//   intPB11.def_property("lastDt",
//                        (Scalar (INT::*)() const) &INT::lastDt,
//                        (void (INT::*)(Scalar)) &INT::lastDt);
//   intPB11.def_property("dtGrowth",
//                        (Scalar (INT::*)() const) &INT::dtGrowth,
//                        (void (INT::*)(Scalar)) &INT::dtGrowth);
//   intPB11.def_property_readonly("dataBase", &INT::dataBase);
//   intPB11.def_property_readonly("physicsPackages", &INT::physicsPackages);
//   intPB11.def_property("rigorousBoundaries",
//                        (bool (INT::*)() const) &INT::rigorousBoundaries,
//                        (void (INT::*)(bool)) &INT::rigorousBoundaries);
//   intPB11.def_property("updateBoundaryFrequency",
//                        (int (INT::*)() const) &INT::updateBoundaryFrequency,
//                        (void (INT::*)(const int)) &INT::updateBoundaryFrequency);
//   intPB11.def_property("verbose",
//                        (bool (INT::*)() const) &INT::verbose,
//                        (void (INT::*)(bool)) &INT::verbose);
//   intPB11.def_property("domainDecompositionIndependent",
//                        (bool (INT::*)() const) &INT::domainDecompositionIndependent,
//                        (void (INT::*)(bool)) &INT::domainDecompositionIndependent);
//   intPB11.def_property("cullGhostNodes",
//                        (bool (INT::*)() const) &INT::cullGhostNodes,
//                        (void (INT::*)(bool)) &INT::cullGhostNodes);
// }

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralGravity, m) {
  using namespace Spheral;

  m.doc() = "Spheral Utilities module: a grab-bag of stuff.";

  //............................................................................
  // KeyTraits
  {
    py::class_<KeyTraits> x(m, "KeyTraits");

    // Attributes
    x.def_readonly_static("numbits", &KeyTraits::numbits);
    x.def_readonly_static("numbits1d", &KeyTraits::numbits1d);
    x.def_readonly_static("zero", &KeyTraits::zero);
    x.def_readonly_static("one", &KeyTraits::one);
    x.def_readonly_static("two", &KeyTraits::two);
    x.def_readonly_static("maxKey1d", &KeyTraits::maxKey1d);
    x.def_readonly_static("maxKey", &KeyTraits::maxKey);
  }

  //............................................................................
  // Timer
  {
    py::class_<Timer> x(m, "Timer");

    // Constructors
    x.def(py::init<>());
    x.def(py::init<std::string>(), "name"_a);
    x.def(py::init<std::string, Timer&>(), "name"_a, "base"_a);

    // Methods
    x.def("setup", &Timer::setup);
    x.def("start", &Timer::start);
    x.def("stop", &Timer::stop);
    x.def("clear", &Timer::clear);
    x.def("getTimeStampWC", &Timer::getTimeStampWC);
    x.def("wc_time", &Timer::wc_time);
    x.def("setup", &Timer::setup);
    x.def("Name", &Timer::Name);
    x.def("Count", &Timer::Count);
    x.def_static("TimerSummary", (void (*)()) &Timer::TimerSummary);
  }

//   //............................................................................
//   // Per dimension bindings.
// #ifdef SPHERAL1D
//   dimensionBindings<Spheral::Dim<1>>(m, "1d");
// #endif

// #ifdef SPHERAL2D
//   dimensionBindings<Spheral::Dim<2>>(m, "2d");
// #endif

// #ifdef SPHERAL3D
//   dimensionBindings<Spheral::Dim<3>>(m, "3d");
// #endif

}
