// Put Python includes first to avoid compile warnings about redefining _POSIX_C_SOURCE
#include "pybind11/pybind11.h"
#include "pybind11/stl_bind.h"
#include "pybind11/operators.h"

#include <vector>
#include <string>

#include "Geometry/Dimension.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/SecondMomentHourglassControl.hh"
#include "Hydro/ThirdMomentHourglassControl.hh"
#include "Hydro/VoronoiHourglassControl.hh"

#include "Pybind11Wraps/DataOutput/PyRestartMethods.hh"
#include "Pybind11Wraps/Physics/PyAbstractPhysics.hh"
#include "Pybind11Wraps/Physics/PyPhysics.hh"
// #include "Pybind11Wraps/Physics/virtualPhysicsBindings.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::PhysicsSpace;

namespace {  // anonymous

// //------------------------------------------------------------------------------
// // Per dimension bindings.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// void dimensionBindings(py::module& m, const std::string suffix) {

//   typedef typename Dimension::Scalar Scalar;
//   typedef typename Dimension::Vector Vector;
//   typedef typename Dimension::Tensor Tensor;
//   typedef typename Dimension::SymTensor SymTensor;
//   typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
//   typedef Spheral::GeomPlane<Dimension> Plane;
//   using Spheral::NodeSpace::NodeList;
//   using Spheral::KernelSpace::TableKernel;

//   //............................................................................
//   // MonaghanGingoldViscosity
//   typedef MonaghanGingoldViscosity<Dimension> MGV;
//   py::class_<MGV,
//              PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, MGV>>> mgPB11(m, ("MonaghanGingoldViscosity" + suffix).c_str());
//   virtualArtificialViscosityBindings<Dimension, MGV>(m, mgPB11);
//   mgPB11
  
//     // Constructors
//     .def(py::init<const Scalar, const Scalar, const bool, const bool>(), "Clinear"_a, "Cquadratic"_a, "linearInExpansion"_a=false, "quadraticInExpansion"_a=false)

//     // Attributes
//     .def_property("linearInExpansion", (bool (MGV::*)() const) &MGV::linearInExpansion, (void (MGV::*)(bool)) &MGV::linearInExpansion)
//     .def_property("quadraticInExpansion", (bool (MGV::*)() const) &MGV::quadraticInExpansion, (void (MGV::*)(bool)) &MGV::quadraticInExpansion)
//     ;

//   //............................................................................
//   // CRKSPHMonaghanGingoldViscosity
//   typedef CRKSPHMonaghanGingoldViscosity<Dimension> CRKMGV;
//   py::class_<CRKMGV, MGV,
//              PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, CRKMGV>>> crkPB11(m, ("CRKSPHMonaghanGingoldViscosity" + suffix).c_str());
//   virtualArtificialViscosityBindings<Dimension, CRKMGV>(m, crkPB11);
//   crkPB11
  
//     // Constructors
//     .def(py::init<const Scalar, const Scalar, const bool, const bool, const Scalar, const Scalar>(),
//          "Clinear"_a=1.0, "Cquadratic"_a=1.0, "linearInExpansion"_a=false, "quadraticInExpansion"_a=false, "etaCritFrac"_a=1.0, "etaFoldFrac"_a=0.2)

//     // Attributes
//     .def_property("etaCritFrac", (Scalar (CRKMGV::*)() const) &CRKMGV::etaCritFrac, (void (CRKMGV::*)(Scalar)) &CRKMGV::etaCritFrac)
//     .def_property("etaFoldFrac", (Scalar (CRKMGV::*)() const) &CRKMGV::etaFoldFrac, (void (CRKMGV::*)(Scalar)) &CRKMGV::etaFoldFrac)
//     ;

//   //............................................................................
//   // MorrisMonaghanReducingViscosity
//   typedef MorrisMonaghanReducingViscosity<Dimension> MMRV;
//   py::class_<MMRV,
//              PhysicsSpace::PyPhysics<Dimension, PhysicsSpace::PyAbstractPhysics<Dimension, MMRV>>> mmrvPB11(m, ("MorrisMonaghanReducingViscosity" + suffix).c_str());
//   // PhysicsSpace::virtualPhysicsBindings<Dimension, MMRV>(m, mmrvPB11);
//   mmrvPB11
  
//     // Constructors
//     .def(py::init<AV&, const Scalar, const Scalar, const Scalar, const Scalar>(),
//          "q"_a, "nhQ"_a=5.0, "nhL"_a=10.0, "aMin"_a=0.1, "aMax"_a=2.0)

//     // Attributes
//     .def_property("nhQ", (Scalar (MMRV::*)() const) &MMRV::nhQ, (void (MMRV::*)(Scalar)) &MMRV::nhQ)
//     .def_property("nhL", (Scalar (MMRV::*)() const) &MMRV::nhL, (void (MMRV::*)(Scalar)) &MMRV::nhL)
//     .def_property("aMin", (Scalar (MMRV::*)() const) &MMRV::aMin, (void (MMRV::*)(Scalar)) &MMRV::aMin)
//     .def_property("aMax", (Scalar (MMRV::*)() const) &MMRV::aMax, (void (MMRV::*)(Scalar)) &MMRV::aMax)
//     .def_property_readonly("DrvAlphaDtQ", &MMRV::DrvAlphaDtQ, py::return_value_policy::reference_internal)
//     .def_property_readonly("DrvAlphaDtL", &MMRV::DrvAlphaDtL, py::return_value_policy::reference_internal)
//     ;

//   //............................................................................
//   // CullenDehnenViscosity
//   typedef CullenDehnenViscosity<Dimension> CDV;
//   py::class_<CDV,
//              PhysicsSpace::PyPhysics<Dimension, PhysicsSpace::PyAbstractPhysics<Dimension, CDV>>> cdvPB11(m, ("CullenDehnenViscosity" + suffix).c_str());
//   cdvPB11
  
//     // Constructors
//     .def(py::init<AV&, const TableKernel<Dimension>&, const Scalar, const Scalar, const Scalar, const Scalar, const Scalar, const Scalar, const bool>(),
//          "q"_a, "W"_a, "alphMax"_a=2.0, "alphMin"_a=0.02, "betaC"_a=0.7, "betaD"_a=0.05, "betaE"_a=1.0, "fKern"_a=0.333333333, "boolHopkins"_a=true)

//     // Attributes
//     .def_property("alphMax", (Scalar (CDV::*)() const) &CDV::alphMax, (void (CDV::*)(Scalar)) &CDV::alphMax)
//     .def_property("alphMin", (Scalar (CDV::*)() const) &CDV::alphMin, (void (CDV::*)(Scalar)) &CDV::alphMin)
//     .def_property("betaE", (Scalar (CDV::*)() const) &CDV::betaE, (void (CDV::*)(Scalar)) &CDV::betaE)
//     .def_property("betaD", (Scalar (CDV::*)() const) &CDV::betaD, (void (CDV::*)(Scalar)) &CDV::betaD)
//     .def_property("betaC", (Scalar (CDV::*)() const) &CDV::betaC, (void (CDV::*)(Scalar)) &CDV::betaC)
//     .def_property("fKern", (Scalar (CDV::*)() const) &CDV::fKern, (void (CDV::*)(Scalar)) &CDV::fKern)
//     .def_property("boolHopkins", (bool (CDV::*)() const) &CDV::boolHopkins, (void (CDV::*)(bool)) &CDV::boolHopkins)
//     .def_property_readonly("PrevDvDt", &CDV::PrevDvDt, py::return_value_policy::reference_internal)
//     .def_property_readonly("PrevDivV", &CDV::PrevDivV, py::return_value_policy::reference_internal)
//     .def_property_readonly("PrevDivV2", &CDV::PrevDivV2, py::return_value_policy::reference_internal)
//     .def_property_readonly("CullAlpha", &CDV::CullAlpha, py::return_value_policy::reference_internal)
//     .def_property_readonly("CullAlpha2", &CDV::CullAlpha2, py::return_value_policy::reference_internal)
//     .def_property_readonly("DalphaDt", &CDV::DalphaDt, py::return_value_policy::reference_internal)
//     .def_property_readonly("alphaLocal", &CDV::alphaLocal, py::return_value_policy::reference_internal)
//     ;

//   //............................................................................
//   // TensorMonaghanGingoldViscosity
//   typedef TensorMonaghanGingoldViscosity<Dimension> TMGV;
//   py::class_<TMGV,
//              PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, TMGV>>> tmgPB11(m, ("TensorMonaghanGingoldViscosity" + suffix).c_str());
//   virtualArtificialViscosityBindings<Dimension, TMGV>(m, tmgPB11);
//   tmgPB11
  
//     // Constructors
//     .def(py::init<const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0)
//     ;

//   //............................................................................
//   // 
//   typedef FiniteVolumeViscosity<Dimension> FVV;
//   py::class_<FVV,
//              PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, FVV>>> fvvPB11(m, ("FiniteVolumeViscosity" + suffix).c_str());
//   virtualArtificialViscosityBindings<Dimension, FVV>(m, fvvPB11);
//   fvvPB11
  
//     // Constructors
//     .def(py::init<const Scalar, const Scalar, const bool>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0, "scalar"_a=false)

//     // Attributes
//     .def_property("scalar", (bool (FVV::*)() const) &FVV::scalar, (void (FVV::*)(bool)) &FVV::scalar)
//     .def_property_readonly("DvDx", &FVV::DvDx, py::return_value_policy::reference_internal)
//     ;

//   //............................................................................
//   // TensorSVPHViscosity
//   typedef TensorSVPHViscosity<Dimension> TSVPHV;
//   py::class_<TSVPHV,
//              PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, TSVPHV>>> tsvphvPB11(m, ("TensorSVPHViscosity" + suffix).c_str());
//   virtualArtificialViscosityBindings<Dimension, TSVPHV>(m, tsvphvPB11);
//   tsvphvPB11
  
//     // Constructors
//     .def(py::init<const Scalar, const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0, "fslice"_a=0.5)

//     // Attributes
//     .def_property("fslice", (Scalar (TSVPHV::*)() const) &TSVPHV::fslice, (void (TSVPHV::*)(Scalar)) &TSVPHV::fslice)
//     .def_property_readonly("shearCorrection", &TSVPHV::shearCorrection, py::return_value_policy::reference_internal)
//     .def_property_readonly("Qface", &TSVPHV::Qface, py::return_value_policy::reference_internal)
//     ;

//   //............................................................................
//   // TensorCRKSPHViscosity
//   typedef TensorCRKSPHViscosity<Dimension> TCRKSPHV;
//   py::class_<TCRKSPHV,
//              PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, TCRKSPHV>>> tcrksphvPB11(m, ("TensorCRKSPHViscosity" + suffix).c_str());
//   virtualArtificialViscosityBindings<Dimension, TCRKSPHV>(m, tcrksphvPB11);
//   tcrksphvPB11
  
//     // Constructors
//     .def(py::init<const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0)
//     ;

//   //............................................................................
//   // VonNeumanViscosity
//   typedef VonNeumanViscosity<Dimension> vNV;
//   py::class_<vNV,
//              PyArtificialViscosity<Dimension, PyAbstractArtificialViscosity<Dimension, vNV>>> vNVPB11(m, ("VonNeumanViscosity" + suffix).c_str());
//   virtualArtificialViscosityBindings<Dimension, vNV>(m, vNVPB11);
//   vNVPB11
  
//     // Constructors
//     .def(py::init<const Scalar, const Scalar>(), "Clinear"_a=1.0, "Cquadratic"_a=1.0)

//     // Attributes
//     .def_property_readonly("viscousEnergy", &vNV::viscousEnergy, py::return_value_policy::reference_internal)
//     ;
// }

} // anonymous

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_PLUGIN(SpheralHydro) {
  using namespace Spheral;
  using namespace Spheral::PhysicsSpace;

  py::module m("SpheralHydro", "Spheral Hydro module.");

  //............................................................................
  // HydroFieldNames
  py::class_<HydroFieldNames>(m, "HydroFieldNames")

    // Static attributes
    .def_readonly_static("mass", &HydroFieldNames::mass)
    .def_readonly_static("position", &HydroFieldNames::position)
    .def_readonly_static("velocity", &HydroFieldNames::velocity)
    .def_readonly_static("H", &HydroFieldNames::H)
    .def_readonly_static("work", &HydroFieldNames::work)
    .def_readonly_static("velocityGradient", &HydroFieldNames::velocityGradient)
    .def_readonly_static("internalVelocityGradient", &HydroFieldNames::internalVelocityGradient)
    .def_readonly_static("massDensity", &HydroFieldNames::massDensity)
    .def_readonly_static("normalization", &HydroFieldNames::normalization)
    .def_readonly_static("specificThermalEnergy", &HydroFieldNames::specificThermalEnergy)
    .def_readonly_static("maxViscousPressure", &HydroFieldNames::maxViscousPressure)
    .def_readonly_static("effectiveViscousPressure", &HydroFieldNames::effectiveViscousPressure)
    .def_readonly_static("viscousWork", &HydroFieldNames::viscousWork)
    .def_readonly_static("XSPHDeltaV", &HydroFieldNames::XSPHDeltaV)
    .def_readonly_static("XSPHWeightSum", &HydroFieldNames::XSPHWeightSum)
    .def_readonly_static("Hsmooth", &HydroFieldNames::Hsmooth)
    .def_readonly_static("massFirstMoment", &HydroFieldNames::massFirstMoment)
    .def_readonly_static("massSecondMoment", &HydroFieldNames::massSecondMoment)
    .def_readonly_static("weightedNeighborSum", &HydroFieldNames::weightedNeighborSum)
    .def_readonly_static("pressure", &HydroFieldNames::pressure)
    .def_readonly_static("temperature", &HydroFieldNames::temperature)
    .def_readonly_static("soundSpeed", &HydroFieldNames::soundSpeed)
    .def_readonly_static("pairAccelerations", &HydroFieldNames::pairAccelerations)
    .def_readonly_static("pairWork", &HydroFieldNames::pairWork)
    .def_readonly_static("omegaGradh", &HydroFieldNames::omegaGradh)
    .def_readonly_static("gamma", &HydroFieldNames::gamma)
    .def_readonly_static("entropy", &HydroFieldNames::entropy)
    .def_readonly_static("PSPHcorrection", &HydroFieldNames::PSPHcorrection)
    .def_readonly_static("numberDensitySum", &HydroFieldNames::numberDensitySum)
    .def_readonly_static("timeStepMask", &HydroFieldNames::timeStepMask)
    .def_readonly_static("m0_CRKSPH", &HydroFieldNames::m0_CRKSPH)
    .def_readonly_static("m1_CRKSPH", &HydroFieldNames::m1_CRKSPH)
    .def_readonly_static("m2_CRKSPH", &HydroFieldNames::m2_CRKSPH)
    .def_readonly_static("m3_CRKSPH", &HydroFieldNames::m3_CRKSPH)
    .def_readonly_static("m4_CRKSPH", &HydroFieldNames::m4_CRKSPH)
    .def_readonly_static("gradM0_CRKSPH", &HydroFieldNames::gradM0_CRKSPH)
    .def_readonly_static("gradM1_CRKSPH", &HydroFieldNames::gradM1_CRKSPH)
    .def_readonly_static("gradM2_CRKSPH", &HydroFieldNames::gradM2_CRKSPH)
    .def_readonly_static("gradM3_CRKSPH", &HydroFieldNames::gradM3_CRKSPH)
    .def_readonly_static("gradM4_CRKSPH", &HydroFieldNames::gradM4_CRKSPH)
    .def_readonly_static("A0_CRKSPH", &HydroFieldNames::A0_CRKSPH)
    .def_readonly_static("A_CRKSPH", &HydroFieldNames::A_CRKSPH)
    .def_readonly_static("B_CRKSPH", &HydroFieldNames::B_CRKSPH)
    .def_readonly_static("C_CRKSPH", &HydroFieldNames::C_CRKSPH)
    .def_readonly_static("gradA_CRKSPH", &HydroFieldNames::gradA_CRKSPH)
    .def_readonly_static("gradB_CRKSPH", &HydroFieldNames::gradB_CRKSPH)
    .def_readonly_static("gradC_CRKSPH", &HydroFieldNames::gradC_CRKSPH)
    .def_readonly_static("surfacePoint", &HydroFieldNames::surfacePoint)
    .def_readonly_static("M_SPHCorrection", &HydroFieldNames::M_SPHCorrection)
    .def_readonly_static("volume", &HydroFieldNames::volume)
    .def_readonly_static("linearMomentum", &HydroFieldNames::linearMomentum)
    .def_readonly_static("totalEnergy", &HydroFieldNames::totalEnergy)
    .def_readonly_static("mesh", &HydroFieldNames::mesh)
    .def_readonly_static("hourglassMask", &HydroFieldNames::hourglassMask)
    .def_readonly_static("faceVelocity", &HydroFieldNames::faceVelocity)
    .def_readonly_static("faceForce", &HydroFieldNames::faceForce)
    .def_readonly_static("faceMass", &HydroFieldNames::faceMass)
    .def_readonly_static("polyvols", &HydroFieldNames::polyvols)
    .def_readonly_static("massDensityGradient", &HydroFieldNames::massDensityGradient)
    .def_readonly_static("ArtificialViscousClMultiplier", &HydroFieldNames::ArtificialViscousClMultiplier)
    .def_readonly_static("ArtificialViscousCqMultiplier", &HydroFieldNames::ArtificialViscousCqMultiplier)
    ;

  //............................................................................
  // Per dimension bindings.
// #ifdef SPHERAL1D
//   dimensionBindings<Spheral::Dim<1>>(m, "1d");
// #endif

// #ifdef SPHERAL2D
//   dimensionBindings<Spheral::Dim<2>>(m, "2d");
//   //............................................................................
//   // MonaghanGingoldViscosityGSRZ
//   py::class_<MonaghanGingoldViscosityGSRZ,
//              PyArtificialViscosity<Dim<2>, PyAbstractArtificialViscosity<Dim<2>, MonaghanGingoldViscosityGSRZ>>> mggsrzPB11(m, "MonaghanGingoldViscosityGSRZ");
//   virtualArtificialViscosityBindings<Dim<2>, MonaghanGingoldViscosityGSRZ>(m, mggsrzPB11);
//   mggsrzPB11
  
//     // Constructors
//     .def(py::init<const Dim<2>::Scalar, const Dim<2>::Scalar, const bool, const bool>(),
//          "Clinear"_a=1.0, "Cquadratic"_a=1.0, "linearInExpansion"_a =false, "quadraticInExpansion"_a=false)
//     ;

// #endif

// #ifdef SPHERAL3D
//   dimensionBindings<Spheral::Dim<3>>(m, "3d");
// #endif

  return m.ptr();
}
