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
#include "Pybind11Wraps/Physics/virtualPhysicsBindings.hh"

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Spheral;
using namespace Spheral::PhysicsSpace;

//------------------------------------------------------------------------------
// Make the module
//------------------------------------------------------------------------------
PYBIND11_MODULE(SpheralHydro, m) {
  using namespace Spheral;
  using namespace Spheral::PhysicsSpace;

  m.doc() = "Spheral Hydro module.";

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
    .def_readonly_static("hydroAcceleration", &HydroFieldNames::hydroAcceleration)
    .def_readonly_static("massDensity", &HydroFieldNames::massDensity)
    .def_readonly_static("normalization", &HydroFieldNames::normalization)
    .def_readonly_static("specificThermalEnergy", &HydroFieldNames::specificThermalEnergy)
    .def_readonly_static("maxViscousPressure", &HydroFieldNames::maxViscousPressure)
    .def_readonly_static("effectiveViscousPressure", &HydroFieldNames::effectiveViscousPressure)
    .def_readonly_static("massDensityCorrection", &HydroFieldNames::massDensityCorrection)
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
    .def_readonly_static("gamma", &HydroFieldNames::gamma)
    .def_readonly_static("entropy", &HydroFieldNames::entropy)
    .def_readonly_static("PSPHcorrection", &HydroFieldNames::PSPHcorrection)
    .def_readonly_static("omegaGradh", &HydroFieldNames::omegaGradh)
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
    .def_readonly_static("gradA_CRKSPH", &HydroFieldNames::gradA0_CRKSPH)
    .def_readonly_static("gradA_CRKSPH", &HydroFieldNames::gradA_CRKSPH)
    .def_readonly_static("gradB_CRKSPH", &HydroFieldNames::gradB_CRKSPH)
    .def_readonly_static("gradC_CRKSPH", &HydroFieldNames::gradC_CRKSPH)
    .def_readonly_static("surfacePoint", &HydroFieldNames::surfacePoint)
    .def_readonly_static("voidPoint", &HydroFieldNames::voidPoint)
    .def_readonly_static("etaVoidPoints", &HydroFieldNames::etaVoidPoints)
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
}
