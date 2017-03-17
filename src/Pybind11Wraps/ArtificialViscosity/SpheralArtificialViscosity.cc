// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosity.hh"
#include "ArtificialViscosity/CRKSPHMonaghanGingoldViscosity.hh"
#include "ArtificialViscosity/MorrisMonaghanReducingViscosity.hh"
#include "ArtificialViscosity/CullenDehnenViscosity.hh"
#include "ArtificialViscosity/TensorMonaghanGingoldViscosity.hh"
#include "ArtificialViscosity/FiniteVolumeViscosity.hh"
#include "ArtificialViscosity/TensorSVPHViscosity.hh"
#include "ArtificialViscosity/TensorCRKSPHViscosity.hh"
#include "ArtificialViscosity/VonNeumanViscosity.hh"
#include "ArtificialViscosity/MonaghanGingoldViscosityGSRZ.hh"

#include "PyAbstractArtificialViscosity.hh"
#include "PyArtificialViscosity.hh"
#include "PyGenericHydro.hh"
#include "PyGenericBodyForce.hh"
#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::ArtificialViscositySpace;

namespace {  // anonymous

//------------------------------------------------------------------------------
// Common virtual methods of ArtificialViscosity objects.
//------------------------------------------------------------------------------
template<typename Dimension, typename Obj, typename PB11Obj>
void virtualArtificialViscosityBindings(py::module& m, PB11Obj& obj) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  obj
    .def("initialize", &Obj::initialize, "dataBase"_a, "state"_a, "derivs"_a, "boundaryBegin"_a, "boundaryEnd"_a, "time"_a, "dt"_a, "W"_a)
    .def("Piij", &Obj::Piij, "xi"_a, "etai"_a, "vi"_a, "rhoi"_a, "csi"_a, "Hi"_a, "xj"_a, "etaj"_a, "vj"_a, "rhoj"_a, "csj"_a, "Hj"_a)
    .def("registerState", &Obj::registerState, "dataBase"_a, "state"_a)
    .def("registerDerivatives", &Obj::registerDerivatives, "dataBase"_a, "state"_a)
    .def("label", &Obj::label)
    .def("dumpState", &Obj::dumpState, "file"_a, "pathName"_a)
    .def("restoreState", &Obj::dumpState, "file"_a, "pathName"_a)
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
  // ArtificialViscosity
  typedef ArtificialViscosity<Dimension> AV;
  py::class_<AV,
             PyAbstractArtificialViscosity<Dimension, AV>> avPB11(m, ("ArtificialViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, Phys>(m, avPB11);
  avPB11
    
    // Constructors
    .def(py::init<>(const Scalar, const Scalar, const CRKSPHSpace::CRKorder), "Clinear"_a, "Cquadratic"_a, "QcorrectionOrder"_a=CRKSPHSpace::CRKOrder::LinearOrder)

    // Methods
    .def("appendBoundary", &Phys::appendBoundary, "boundary"_a)
    .def("prependBoundary", &Phys::prependBoundary, "boundary"_a)
    .def("clearBoundaries", &Phys::clearBoundaries)
    .def("haveBoundary", &Phys::haveBoundary, "boundary"_a)
    .def("boundaryConditions", &Phys::boundaryConditions)

    // Attributes
    .def_property("Cl", (Scalar (AV::*)() const) &AV::Cl, (void (AV::*)(Scalar)) &AV::Cl)
    ;
}

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralArtificialViscosity) {
  using namespace Spheral;
  using namespace Spheral::ArtificialViscositySpace;

  py::module m("SpheralArtificialViscosity", "Spheral ArtificialViscosity module.");

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

  return m.ptr();
}
