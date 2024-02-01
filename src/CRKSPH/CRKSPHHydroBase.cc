//---------------------------------Spheral++----------------------------------//
// Hydro -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "RK/ReproducingKernel.hh"
#include "RK/RKFieldNames.hh"
#include "CRKSPH/computeCRKSPHSumMassDensity.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Hydro/entropyWeightingFunction.hh"
#include "Strength/SolidFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/updateStateFields.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
#include "Hydro/EntropyPolicy.hh"
#include "RK/ContinuityVolumePolicy.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/newtonRaphson.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/computeShepardsInterpolation.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"

#include "CRKSPHHydroBase.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

using std::vector;
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
template<typename Dimension>
CRKSPHHydroBase<Dimension>::
CRKSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                DataBase<Dimension>& dataBase,
                ArtificialViscosity<Dimension>& Q,
                const RKOrder order,
                const double filter,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool evolveTotalEnergy,
                const bool XSPH,
                const MassDensityType densityUpdate,
                const HEvolutionType HUpdate,
                const double epsTensile,
                const double nTensile):
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mOrder(order),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mfilter(filter),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mSpecificThermalEnergy0(FieldStorageType::CopyFields),
  mEntropy(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
  mMaxViscousPressure(FieldStorageType::CopyFields),
  mEffViscousPressure(FieldStorageType::CopyFields),
  mViscousWork(FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldStorageType::CopyFields),
  mMassSecondMoment(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mPairAccelerations(),
  mRestart(registerWithRestart(*this)) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  mEntropy = dataBase.newFluidFieldList(0.0, HydroFieldNames::entropy);
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H);
  mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
  mEffViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::effectiveViscousPressure);
  mViscousWork = dataBase.newFluidFieldList(0.0, HydroFieldNames::viscousWork);
  mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::H);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mPairAccelerations.clear();
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHHydroBase<Dimension>::
~CRKSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
initializeProblemStartupDependencies(DataBase<Dimension>& dataBase,
                                     State<Dimension>& state,
                                     StateDerivatives<Dimension>& derivs) {

  // Initialize the pressure, sound speed, and entropy.
  updateStateFields<Dimension>(HydroFieldNames::pressure, state, derivs);
  updateStateFields<Dimension>(HydroFieldNames::soundSpeed, state, derivs);
  updateStateFields<Dimension>(HydroFieldNames::entropy, state, derivs);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // Create the local storage for time step mask, pressure, sound speed, and correction fields.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeFluidFieldList(mEntropy,    0.0,                   HydroFieldNames::entropy, false);
  dataBase.resizeFluidFieldList(mPressure,   0.0,                   HydroFieldNames::pressure, false);
  dataBase.resizeFluidFieldList(mSoundSpeed, 0.0,                   HydroFieldNames::soundSpeed, false);

  // We have to choose either compatible or total energy evolution.
  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "CRKSPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // If we're using the compatibile energy discretization, prepare to maintain a copy
  // of the thermal energy.
  dataBase.resizeFluidFieldList(mSpecificThermalEnergy0, 0.0);
  auto nodeListi = 0u;
  if (mCompatibleEnergyEvolution) {
    for (auto itr = dataBase.fluidNodeListBegin(); itr < dataBase.fluidNodeListEnd(); ++itr, ++nodeListi) {
      *mSpecificThermalEnergy0[nodeListi] = (*itr)->specificThermalEnergy();
      (*mSpecificThermalEnergy0[nodeListi]).name(HydroFieldNames::specificThermalEnergy + "0");
    }
  }

  // Now register away.
  // Mass.
  auto mass = dataBase.fluidMass();
  state.enroll(mass);

  // Volume.
  // Note: we depend on RKCorrections having already registered the volume, but it was registered
  // without an update policy.
  auto vol = state.fields(HydroFieldNames::volume, 0.0);
  CHECK(vol.size() == dataBase.numFluidNodeLists());
  state.enroll(vol, make_policy<ContinuityVolumePolicy<Dimension>>());

  // We need to build up CompositeFieldListPolicies for the mass density and H fields
  // in order to enforce NodeList dependent limits.
  auto massDensity = dataBase.fluidMassDensity();
  auto Hfield = dataBase.fluidHfield();
  nodeListi = 0u;
  for (auto itr = dataBase.fluidNodeListBegin();
       itr < dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    state.enroll(*massDensity[nodeListi], make_policy<IncrementBoundedState<Dimension, Scalar>>((*itr)->rhoMin(),
                                                                                                (*itr)->rhoMax()));
    const auto hmaxInv = 1.0/(*itr)->hmax();
    const auto hminInv = 1.0/(*itr)->hmin();
    switch (this->HEvolution()) {
      case HEvolutionType::IntegrateH:
        state.enroll(*Hfield[nodeListi], make_policy<IncrementBoundedState<Dimension, SymTensor, Scalar>>(hmaxInv, hminInv));
        break;

      case HEvolutionType::IdealH:
        state.enroll(*Hfield[nodeListi], make_policy<ReplaceBoundedState<Dimension, SymTensor, Scalar>>(hmaxInv, hminInv));
        break;

       default:
         VERIFY2(false, "SPH ERROR: Unknown Hevolution option ");
    }
  }

  // Register the position update, which depends on whether we're using XSPH or not.
  auto position = dataBase.fluidPosition();
  state.enroll(position, make_policy<IncrementState<Dimension, Vector>>());

  // Register the entropy.
  state.enroll(mEntropy, make_policy<EntropyPolicy<Dimension>>());

  // Are we using the compatible energy evolution scheme?
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  auto velocity = dataBase.fluidVelocity();
  if (compatibleEnergyEvolution()) {
    // The compatible energy update.
    state.enroll(specificThermalEnergy, make_policy<SpecificThermalEnergyPolicy<Dimension>>(dataBase));
    state.enroll(velocity, make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::specificThermalEnergy},
                                                                          true));
    state.enroll(mSpecificThermalEnergy0);

  } else if (mEvolveTotalEnergy) {
    // If we're doing total energy, we register the specific energy to advance with the
    // total energy policy.
    state.enroll(specificThermalEnergy, make_policy<SpecificFromTotalThermalEnergyPolicy<Dimension>>());
    state.enroll(velocity, make_policy<IncrementState<Dimension, Vector>>({HydroFieldNames::specificThermalEnergy},
                                                                          true));

  } else {
    // Otherwise we're just time-evolving the specific energy.
    state.enroll(specificThermalEnergy, make_policy<IncrementState<Dimension, Scalar>>());
    state.enroll(velocity, make_policy<IncrementState<Dimension, Vector>>(true));
  }

  // Register the time step mask, initialized to 1 so that everything defaults to being
  // checked.
  state.enroll(mTimeStepMask);

  // Compute and register the pressure and sound speed.
  state.enroll(mPressure, make_policy<PressurePolicy<Dimension>>());
  state.enroll(mSoundSpeed, make_policy<SoundSpeedPolicy<Dimension>>());
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  const auto DxDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const auto DvDtName = HydroFieldNames::hydroAcceleration;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, SymTensor>::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  dataBase.resizeFluidFieldList(mEffViscousPressure, 0.0, HydroFieldNames::effectiveViscousPressure, false);
  dataBase.resizeFluidFieldList(mViscousWork, 0.0, HydroFieldNames::viscousWork, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);

  derivs.enroll(mHideal);
  derivs.enroll(mMaxViscousPressure);
  derivs.enroll(mEffViscousPressure);
  derivs.enroll(mViscousWork);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mXSPHDeltaV);

  // These two (the position and velocity updates) may be registered
  // by other physics packages as well, so we need to be careful
  // not to duplicate if so.
  if (not derivs.registered(mDxDt)) derivs.enroll(mDxDt);
  if (not derivs.registered(mDvDt)) derivs.enroll(mDvDt);

  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDHDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enrollAny(HydroFieldNames::pairAccelerations, mPairAccelerations);
}

//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  // Note: we depend on RKCorrections having already updated the volume!
  if (mDensityUpdate == MassDensityType::RigorousSumDensity or
      mDensityUpdate == MassDensityType::VoronoiCellDensity) {
    auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const auto& WR = state.template getAny<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(mOrder));
    const auto& W = WR.kernel();
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto  vol = state.fields(HydroFieldNames::volume, 0.0);
    if (densityUpdate() == MassDensityType::RigorousSumDensity) {
      computeCRKSPHSumMassDensity(connectivityMap, W, position, mass, vol, H, massDensity);
    } else {
      massDensity.assignFields(mass/vol);
    }
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr != this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr != this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro before evaluating derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  // Initialize the artificial viscosity
  const auto& WR = state.template getAny<ReproducingKernel<Dimension>>(RKFieldNames::reproducingKernel(mOrder));
  auto&       Q = this->artificialViscosity();
  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               WR.kernel());
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    auto DepsDt = derivs.fields(IncrementState<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) {
      (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
      (*boundaryItr)->applyFieldListGhostBoundary(DepsDt);
    }
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {

  // Apply boundary conditions to the basic fluid state Fields.
  // volume, mass, and massDensity handled by RKCorrections
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto entropy = state.fields(HydroFieldNames::entropy, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) {
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
  }

  for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(entropy);
    if (compatibleEnergyEvolution()) (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  // Enforce boundary conditions on the fluid state Fields.
  // volume, mass, and massDensity handled by RKCorrections
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  auto entropy = state.fields(HydroFieldNames::entropy, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) {
    specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);
  }

  for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->enforceFieldListBoundary(entropy);
    if (compatibleEnergyEvolution()) (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
  }
}

//------------------------------------------------------------------------------
// Return the RK orders we want to use
//------------------------------------------------------------------------------
template<typename Dimension>
std::set<RKOrder>
CRKSPHHydroBase<Dimension>::
requireReproducingKernels() const {
  return std::set<RKOrder>({RKOrder::ZerothOrder, mOrder});
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.write(mEntropy, pathName + "/entropy");
  file.write(mHideal, pathName + "/Hideal");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.write(mEffViscousPressure, pathName + "/effViscousPressure");
  file.write(mViscousWork, pathName + "/viscousWork");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.read(mEntropy, pathName + "/entropy");
  file.read(mHideal, pathName + "/Hideal");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.read(mEffViscousPressure, pathName + "/effViscousPressure");
  file.read(mViscousWork, pathName + "/viscousWork");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
}

}

