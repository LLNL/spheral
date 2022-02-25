//---------------------------------Spheral++----------------------------------//
// SphericalSPHHydroBase -- An SPH/ASPH hydrodynamic package for Spheral++,
//                          specialized for 1D Spherical (r) geometry.
//
// Note this version is currently abusing our ordinary 1D geometric types,
// implicitly mapping x->r.
//
// Created by JMO, Tue Dec 22 10:04:21 PST 2020
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "computeSPHSumMassDensity.hh"
#include "correctSPHSumMassDensity.hh"
#include "computeSumVoronoiCellMassDensity.hh"
#include "computeSPHOmegaGradhCorrection.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "Hydro/VolumePolicy.hh"
#include "Hydro/VoronoiMassDensityPolicy.hh"
#include "Hydro/SumVoronoiMassDensityPolicy.hh"
//#include "Hydro/SphericalNonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "Mesh/MeshPolicy.hh"
#include "Mesh/generateMesh.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Mesh/Mesh.hh"
#include "CRKSPH/volumeSpacing.hh"

#include "SphericalSPHHydroBase.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
using std::vector;
using std::map;
using std::string;
using std::pair;
using std::make_pair;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
SphericalSPHHydroBase::
SphericalSPHHydroBase(const SmoothingScaleBase<Dim<1>>& smoothingScaleMethod,
                      DataBase<Dimension>& dataBase,
                      ArtificialViscosity<Dim<1>>& Q,
                      const SphericalKernel& W,
                      const SphericalKernel& WPi,
                      const double filter,
                      const double cfl,
                      const bool useVelocityMagnitudeForDt,
                      const bool compatibleEnergyEvolution,
                      const bool evolveTotalEnergy,
                      const bool gradhCorrection,
                      const bool XSPH,
                      const bool correctVelocityGradient,
                      const bool sumMassDensityOverAllNodeLists,
                      const MassDensityType densityUpdate,
                      const HEvolutionType HUpdate,
                      const double epsTensile,
                      const double nTensile,
                      const Vector& xmin,
                      const Vector& xmax):
  SPHHydroBase<Dim<1>>(smoothingScaleMethod,
                       dataBase,
                       Q,
                       W.baseKernel1d(),
                       WPi.baseKernel1d(),
                       filter,
                       cfl,
                       useVelocityMagnitudeForDt,
                       compatibleEnergyEvolution,
                       evolveTotalEnergy,
                       gradhCorrection,
                       XSPH,
                       correctVelocityGradient,
                       sumMassDensityOverAllNodeLists,
                       densityUpdate,
                       HUpdate,
                       epsTensile,
                       nTensile,
                       xmin,
                       xmax),
  mKernel(W),
  mPiKernel(WPi) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
SphericalSPHHydroBase::
~SphericalSPHHydroBase() {
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
registerState(DataBase<Dim<1>>& dataBase,
              State<Dim<1>>& state) {

  // typedef State<Dimension>::PolicyPointer PolicyPointer;

  // The base class does most of it.
  SPHHydroBase<Dim<1>>::registerState(dataBase, state);

  // // Are we using the compatible energy evolution scheme?
  // // If so we need to override the ordinary energy registration with a specialized version.
  // if (mCompatibleEnergyEvolution) {
  //   FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  //   PolicyPointer thermalEnergyPolicy(new SphericalNonSymmetricSpecificThermalEnergyPolicy(dataBase));
  //   state.enroll(specificThermalEnergy, thermalEnergyPolicy);

  //   // Get the policy for the position, and add the specific energy as a dependency.
  //   PolicyPointer positionPolicy = state.policy(state.buildFieldKey(HydroFieldNames::position, UpdatePolicyBase<Dimension>::wildcard()));
  //   positionPolicy->addDependency(HydroFieldNames::specificThermalEnergy);
  // }
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // We have to override the density summation to use the correct kernel
  switch(densityUpdate()) {

  case MassDensityType::IntegrateDensity:
    break;

  case MassDensityType::RigorousSumDensity:
  case MassDensityType::CorrectedSumDensity:
    {
      const auto& connectivityMap = dataBase.connectivityMap();
      const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
      const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
      const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
      auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      computeSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      if (densityUpdate() == MassDensityType::CorrectedSumDensity) {
        correctSPHSumMassDensity(connectivityMap, this->kernel(), mSumMassDensityOverAllNodeLists, position, mass, H, massDensity);
        for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
        for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
      }
    }
    break;

  case MassDensityType::SumDensity:
    {
      auto       massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
      const auto massDensitySum = derivs.fields(ReplaceFieldList<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                HydroFieldNames::massDensity, 0.0);
      massDensity.assignFields(massDensitySum);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
      for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
    }
    break;

  default:
    VERIFY2(false, "Unsupported mass density definition for Spherical SPH");
    break;
  }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
evaluateDerivatives(const Dim<1>::Scalar time,
                    const Dim<1>::Scalar dt,
                    const DataBase<Dim<1>>& dataBase,
                    const State<Dim<1>>& state,
                    StateDerivatives<Dim<1>>& derivs) const {
  this->evaluateDerivativesImpl(time,
                                dt,
                                dataBase,
                                state,
                                derivs,
                                mKernel,
                                mPiKernel,
                                mKernel.baseKernel1d());
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
applyGhostBoundaries(State<Dim<1>>& state,
                     StateDerivatives<Dim<1>>& derivs) {

  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) /= dA;
    }
  }

  // Apply ordinary SPH BCs.
  SPHHydroBase<Dim<1>>::applyGhostBoundaries(state, derivs);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Scale back to mass.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
void
SphericalSPHHydroBase::
enforceBoundaries(State<Dim<1>>& state,
                  StateDerivatives<Dim<1>>& derivs) {

  // Convert the mass to mass/length^2 before BCs are applied.
  auto       mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto pos = state.fields(HydroFieldNames::position, Vector::zero);
  const auto numNodeLists = mass.numFields();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) /= dA;
    }
  }

  // Apply ordinary SPH BCs.
  SPHHydroBase<Dim<1>>::enforceBoundaries(state, derivs);

  // Scale back to mass.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto n = mass[nodeListi]->numElements();
    for (auto i = 0u; i < n; ++i) {
      const auto dA = FastMath::pow2(pos(nodeListi, i).x());
      CHECK(dA > 0.0);
      mass(nodeListi, i) *= dA;
    }
  }
}

//------------------------------------------------------------------------------
// Access the main kernel used for (A)SPH field estimates.
//------------------------------------------------------------------------------
const SphericalKernel&
SphericalSPHHydroBase::
kernel() const {
  return mKernel;
}

//------------------------------------------------------------------------------
// Access the kernel used for artificial viscosity gradients.
//------------------------------------------------------------------------------
const SphericalKernel&
SphericalSPHHydroBase::
PiKernel() const {
  return mPiKernel;
}

}
