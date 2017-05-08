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
    .def("Piij", &Obj::Piij, "nodeListi"_a, "i"_a, "nodeListj"_a, "j"_a, "xi"_a, "etai"_a, "vi"_a, "rhoi"_a, "csi"_a, "Hi"_a, "xj"_a, "etaj"_a, "vj"_a, "rhoj"_a, "csj"_a, "Hj"_a)
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
  using Spheral::ArtificialViscositySpace::MonaghanGingoldViscosity;

  //............................................................................
  // ArtificialViscosity
  typedef ArtificialViscosity<Dimension> AV;
  py::class_<AV,
             PyAbstractArtificialViscosity<Dimension, AV>> avPB11(m, ("ArtificialViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, AV>(m, avPB11);
  avPB11
    
    // Constructors
    .def(py::init<const Scalar, const Scalar, const CRKSPHSpace::CRKOrder>(), "Clinear"_a, "Cquadratic"_a, "QcorrectionOrder"_a=CRKSPHSpace::CRKOrder::LinearOrder)

    // Methods
    .def("curlVelocityMagnitude", &AV::curlVelocityMagnitude, "DvDx"_a)

    // Attributes
    .def_property("Cl", (Scalar (AV::*)() const) &AV::Cl, (void (AV::*)(Scalar)) &AV::Cl)
    .def_property("Cq", (Scalar (AV::*)() const) &AV::Cq, (void (AV::*)(Scalar)) &AV::Cq)
    .def_property("QcorrectionOrder", (CRKSPHSpace::CRKOrder (AV::*)() const) &AV::QcorrectionOrder, (void (AV::*)(CRKSPHSpace::CRKOrder)) &AV::QcorrectionOrder)
    .def_property("balsaraShearCorrection", (bool (AV::*)() const) &AV::balsaraShearCorrection, (void (AV::*)(bool)) &AV::balsaraShearCorrection)
    .def_property("limiter", (bool (AV::*)() const) &AV::limiter, (void (AV::*)(bool)) &AV::limiter)
    .def_property("epsilon2", (Scalar (AV::*)() const) &AV::epsilon2, (void (AV::*)(Scalar)) &AV::epsilon2)
    .def_property("negligibleSoundSpeed", (Scalar (AV::*)() const) &AV::negligibleSoundSpeed, (void (AV::*)(Scalar)) &AV::negligibleSoundSpeed)
    .def_property("csMultiplier", (Scalar (AV::*)() const) &AV::csMultiplier, (void (AV::*)(Scalar)) &AV::csMultiplier)
    .def_property("energyMultiplier", (Scalar (AV::*)() const) &AV::energyMultiplier, (void (AV::*)(Scalar)) &AV::energyMultiplier)
    ;

  //............................................................................
  // MonaghanGingoldViscosity
  typedef MonaghanGingoldViscosity<Dimension> MGV;
  py::class_<MGV,
             PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, MGV>>> mgPB11(m, ("MonaghanGingoldViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, MGV>(m, mgPB11);
  mgPB11
  
    // Constructors
    .def(py::init<const Scalar, const Scalar, const bool, const bool>(), "Clinear"_a, "Cquadratic"_a, "linearInExpansion"_a=false, "quadraticInExpansion"_a=false)

    // Attributes
    .def_property("linearInExpansion", (bool (MGV::*)() const) &MGV::linearInExpansion, (void (MGV::*)(bool)) &MGV::linearInExpansion)
    .def_property("quadraticInExpansion", (bool (MGV::*)() const) &MGV::quadraticInExpansion, (void (MGV::*)(bool)) &MGV::quadraticInExpansion)
    ;

  //............................................................................
  // CRKSPHMonaghanGingoldViscosity
  typedef CRKSPHMonaghanGingoldViscosity<Dimension> CRKMGV;
  py::class_<CRKMGV, MGV,
             PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, CRKMGV>>> crkPB11(m, ("CRKSPHMonaghanGingoldViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, CRKMGV>(m, crkPB11);
  crkPB11
  
    // Constructors
    .def(py::init<const Scalar, const Scalar, const bool, const bool, const Scalar, const Scalar>(),
         "Clinear"_a=1.0, "Cquadratic"_a=1.0, "linearInExpansion"_a=false, "quadraticInExpansion"_a=false, "etaCritFrac"_a=1.0, "etaFoldFrac"_a=0.2)

    // Attributes
    .def_property("etaCritFrac", (Scalar (CRKMGV::*)() const) &CRKMGV::etaCritFrac, (void (CRKMGV::*)(Scalar)) &CRKMGV::etaCritFrac)
    .def_property("etaFoldFrac", (Scalar (CRKMGV::*)() const) &CRKMGV::etaFoldFrac, (void (CRKMGV::*)(Scalar)) &CRKMGV::etaFoldFrac)
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
