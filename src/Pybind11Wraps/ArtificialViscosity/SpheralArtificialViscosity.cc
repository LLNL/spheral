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
#include "Pybind11Wraps/Physics/PyAbstractPhysics.hh"
#include "Pybind11Wraps/Physics/PyPhysics.hh"
// #include "Pybind11Wraps/Physics/virtualPhysicsBindings.hh"

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

  //............................................................................
  // MorrisMonaghanReducingViscosity
  typedef MorrisMonaghanReducingViscosity<Dimension> MMRV;
  py::class_<MMRV,
             PhysicsSpace::PyPhysics<Dimension, PhysicsSpace::PyAbstractPhysics<Dimension, MMRV>>> mmrvPB11(m, ("MorrisMonaghanReducingViscosity" + suffix).c_str());
  // PhysicsSpace::virtualPhysicsBindings<Dimension, MMRV>(m, mmrvPB11);
  mmrvPB11
  
    // Constructors
    .def(py::init<AV&, const Scalar, const Scalar, const Scalar, const Scalar>(),
         "q"_a, "nhQ"_a=5.0, "nhL"_a=10.0, "aMin"_a=0.1, "aMax"_a=2.0)

    // Attributes
    .def_property("nhQ", (Scalar (MMRV::*)() const) &MMRV::nhQ, (void (MMRV::*)(Scalar)) &MMRV::nhQ)
    .def_property("nhL", (Scalar (MMRV::*)() const) &MMRV::nhL, (void (MMRV::*)(Scalar)) &MMRV::nhL)
    .def_property("aMin", (Scalar (MMRV::*)() const) &MMRV::aMin, (void (MMRV::*)(Scalar)) &MMRV::aMin)
    .def_property("aMax", (Scalar (MMRV::*)() const) &MMRV::aMax, (void (MMRV::*)(Scalar)) &MMRV::aMax)
    .def_property_readonly("DrvAlphaDtQ", &MMRV::DrvAlphaDtQ, py::return_value_policy::reference_internal)
    .def_property_readonly("DrvAlphaDtL", &MMRV::DrvAlphaDtL, py::return_value_policy::reference_internal)
    ;

  //............................................................................
  // CullenDehnenViscosity
  typedef CullenDehnenViscosity<Dimension> CDV;
  py::class_<CDV,
             PhysicsSpace::PyPhysics<Dimension, PhysicsSpace::PyAbstractPhysics<Dimension, CDV>>> cdvPB11(m, ("CullenDehnenViscosity" + suffix).c_str());
  cdvPB11
  
    // Constructors
    .def(py::init<AV&, const TableKernel<Dimension>&, const Scalar, const Scalar, const Scalar, const Scalar, const Scalar, const Scalar, const bool>(),
         "q"_a, "W"_a, "alphMax"_a=2.0, "alphMin"_a=0.02, "betaC"_a=0.7, "betaD"_a=0.05, "betaE"_a=1.0, "fKern"_a=0.333333333, "boolHopkins"_a=true)

    // Attributes
    .def_property("alphMax", (Scalar (CDV::*)() const) &CDV::alphMax, (void (CDV::*)(Scalar)) &CDV::alphMax)
    .def_property("alphMin", (Scalar (CDV::*)() const) &CDV::alphMin, (void (CDV::*)(Scalar)) &CDV::alphMin)
    .def_property("betaE", (Scalar (CDV::*)() const) &CDV::betaE, (void (CDV::*)(Scalar)) &CDV::betaE)
    .def_property("betaD", (Scalar (CDV::*)() const) &CDV::betaD, (void (CDV::*)(Scalar)) &CDV::betaD)
    .def_property("betaC", (Scalar (CDV::*)() const) &CDV::betaC, (void (CDV::*)(Scalar)) &CDV::betaC)
    .def_property("fKern", (Scalar (CDV::*)() const) &CDV::fKern, (void (CDV::*)(Scalar)) &CDV::fKern)
    .def_property("boolHopkins", (bool (CDV::*)() const) &CDV::boolHopkins, (void (CDV::*)(bool)) &CDV::boolHopkins)
    .def_property_readonly("PrevDvDt", &CDV::PrevDvDt, py::return_value_policy::reference_internal)
    .def_property_readonly("PrevDivV", &CDV::PrevDivV, py::return_value_policy::reference_internal)
    .def_property_readonly("PrevDivV2", &CDV::PrevDivV2, py::return_value_policy::reference_internal)
    .def_property_readonly("CullAlpha", &CDV::CullAlpha, py::return_value_policy::reference_internal)
    .def_property_readonly("CullAlpha2", &CDV::CullAlpha2, py::return_value_policy::reference_internal)
    .def_property_readonly("DalphaDt", &CDV::DalphaDt, py::return_value_policy::reference_internal)
    .def_property_readonly("alphaLocal", &CDV::alphaLocal, py::return_value_policy::reference_internal)
    ;

  //............................................................................
  // TensorMonaghanGingoldViscosity
  typedef TensorMonaghanGingoldViscosity<Dimension> TMGV;
  py::class_<TMGV,
             PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, TMGV>>> tmgPB11(m, ("TensorMonaghanGingoldViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, TMGV>(m, tmgPB11);
  tmgPB11
  
    // Constructors
    .def(py::init<const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0)
    ;

  //............................................................................
  // 
  typedef FiniteVolumeViscosity<Dimension> FVV;
  py::class_<FVV,
             PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, FVV>>> fvvPB11(m, ("FiniteVolumeViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, FVV>(m, fvvPB11);
  fvvPB11
  
    // Constructors
    .def(py::init<const Scalar, const Scalar, const bool>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0, "scalar"_a=false)

    // Attributes
    .def_property("scalar", (bool (FVV::*)() const) &FVV::scalar, (void (FVV::*)(bool)) &FVV::scalar)
    .def_property_readonly("DvDx", &FVV::DvDx, py::return_value_policy::reference_internal)
    ;

  //............................................................................
  // TensorSVPHViscosity
  typedef TensorSVPHViscosity<Dimension> TSVPHV;
  py::class_<TSVPHV,
             PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, TSVPHV>>> tsvphvPB11(m, ("TensorSVPHViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, TSVPHV>(m, tsvphvPB11);
  tsvphvPB11
  
    // Constructors
    .def(py::init<const Scalar, const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0, "fslice"_a=0.5)

    // Attributes
    .def_property("fslice", (Scalar (TSVPHV::*)() const) &TSVPHV::fslice, (void (TSVPHV::*)(Scalar)) &TSVPHV::fslice)
    .def_property_readonly("shearCorrection", &TSVPHV::shearCorrection, py::return_value_policy::reference_internal)
    .def_property_readonly("Qface", &TSVPHV::Qface, py::return_value_policy::reference_internal)
    ;

  //............................................................................
  // TensorCRKSPHViscosity
  typedef TensorCRKSPHViscosity<Dimension> TCRKSPHV;
  py::class_<TCRKSPHV,
             PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, TCRKSPHV>>> tcrksphvPB11(m, ("TensorCRKSPHViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, TCRKSPHV>(m, tcrksphvPB11);
  tcrksphvPB11
  
    // Constructors
    .def(py::init<const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0)
    ;

  //............................................................................
  // VonNeumanViscosity
  typedef VonNeumanViscosity<Dimension> vNV;
  py::class_<vNV,
             PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, vNV>>> vNVPB11(m, ("VonNeumanViscosity" + suffix).c_str());
  virtualArtificialViscosityBindings<Dimension, vNV>(m, vNVPB11);
  vNVPB11
  
    // Constructors
    .def(py::init<const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0)

    // Attributes
    .def_property_readonly("viscousEnergy", &vNV::viscousEnergy, py::return_value_policy::reference_internal)
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

#ifdef SPHERAL2D
  dimensionBindings<Spheral::Dim<2>>(m, "2d");
  //............................................................................
  // MonaghanGingoldViscosityGSRZ
  py::class_<MonaghanGingoldViscosityGSRZ,
             PyArtificialViscosity<Dim<2>, PyAbstractArtificialViscosity<Dim<2>, MonaghanGingoldViscosityGSRZ>>> mggsrzPB11(m, "MonaghanGingoldViscosityGSRZ");
  virtualArtificialViscosityBindings<Dim<2>, MonaghanGingoldViscosityGSRZ>(m, mggsrzPB11);
  mggsrzPB11
  
    // Constructors
    .def(py::init<const Dim<2>::Scalar, const Dim<2>::Scalar, const bool, const bool>(),
         "Clinear"_a=1.0, "Cquadratic"_a=1.0, "linearInExpansion"_a =false, "quadraticInExpansion"_a=false)
    ;

#endif

#ifdef SPHERAL3D
  dimensionBindings<Spheral::Dim<3>>(m, "3d");
#endif

  return m.ptr();
}
