//---------------------------------Spheral++----------------------------------//
// Hydro -- The CRKSPH/ACRKSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "CRKSPHHydroBase.hh"
#include "CRKSPHUtilities.hh"
#include "computeHullVolumes.hh"
#include "computeNeighborHull.hh"
#include "computeCRKSPHSumMassDensity.hh"
#include "computeHullSumMassDensity.hh"
#include "computeCRKSPHCorrections.hh"
#include "computeCRKSPHIntegral.hh"
#include "computeHVolumes.hh"
#include "centerOfMass.hh"
#include "computeVoronoiCentroids.hh"
#include "volumeSpacing.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"
#include "CRKSPHSpecificThermalEnergyPolicy.hh"
#include "HVolumePolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/NonSymmetricSpecificThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
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
#include "FileIO/FileIO.hh"

#include "SPH/computeSPHSumMassDensity.hh"
#include "gradientCRKSPH.hh"
#include "Geometry/innerProduct.hh"
#include "Geometry/outerProduct.hh"

namespace Spheral {
namespace CRKSPHSpace {

using namespace std;
using PhysicsSpace::Physics;
using PhysicsSpace::GenericHydro;
using NodeSpace::SmoothingScaleBase;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FileIOSpace::FileIO;
using ArtificialViscositySpace::ArtificialViscosity;
using KernelSpace::TableKernel;
using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using NeighborSpace::ConnectivityMap;
using Geometry::innerProduct;
using Geometry::outerProduct;

using PhysicsSpace::MassDensityType;
using PhysicsSpace::HEvolutionType;

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
CRKSPHHydroBase<Dimension>::
CRKSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                const TableKernel<Dimension>& W,
                const TableKernel<Dimension>& WPi,
                ArtificialViscosity<Dimension>& Q,
                const double filter,
                const double cfl,
                const bool useVelocityMagnitudeForDt,
                const bool compatibleEnergyEvolution,
                const bool XSPH,
                const MassDensityType densityUpdate,
                const HEvolutionType HUpdate,
                const double epsTensile,
                const double nTensile):
  GenericHydro<Dimension>(W, WPi, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mXSPH(XSPH),
  mfilter(filter),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mTimeStepMask(FieldSpace::Copy),
  mPressure(FieldSpace::Copy),
  mSoundSpeed(FieldSpace::Copy),
  mVolume(FieldSpace::Copy),
  mSpecificThermalEnergy0(FieldSpace::Copy),
  mHideal(FieldSpace::Copy),
  mMaxViscousPressure(FieldSpace::Copy),
  mEffViscousPressure(FieldSpace::Copy),
  mViscousWork(FieldSpace::Copy),
  mWeightedNeighborSum(FieldSpace::Copy),
  mMassSecondMoment(FieldSpace::Copy),
  mXSPHDeltaV(FieldSpace::Copy),
  mDxDt(FieldSpace::Copy),
  mDvDt(FieldSpace::Copy),
  mDmassDensityDt(FieldSpace::Copy),
  mDspecificThermalEnergyDt(FieldSpace::Copy),
  mDHDt(FieldSpace::Copy),
  mDvDx(FieldSpace::Copy),
  mInternalDvDx(FieldSpace::Copy),
  mDmassDensityDx(FieldSpace::Copy),
  mPairAccelerations(FieldSpace::Copy),
  mA(FieldSpace::Copy),
  mB(FieldSpace::Copy),
  mGradA(FieldSpace::Copy),
  mGradB(FieldSpace::Copy),
  mSurfNorm(FieldSpace::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
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
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
  mEffViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::effectiveViscousPressure);
  mViscousWork = dataBase.newFluidFieldList(0.0, HydroFieldNames::viscousWork);
  mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mDmassDensityDx = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::massDensityGradient);
  mPairAccelerations = dataBase.newFluidFieldList(vector<Vector>(), HydroFieldNames::pairAccelerations);

  mA = dataBase.newFluidFieldList(0.0,              HydroFieldNames::A_CRKSPH);
  mB = dataBase.newFluidFieldList(Vector::zero,     HydroFieldNames::B_CRKSPH);
  mGradA = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::gradA_CRKSPH);
  mGradB = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::gradB_CRKSPH);
  // mSurfNorm = dataBase.newFluidFieldList(Vector::zero, "Surface Normal");

  // // Compute the volumes.
  // const TableKernel<Dimension>& W = this->kernel();
  // const FieldList<Dimension, SymTensor> H = dataBase.fluidHfield();
  // computeHVolumes(W.kernelExtent(), H, mVolume);

  // // We need the boundary information on volume to initialize the CRKSPH corrections.
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->applyFieldListGhostBoundary(mVolume);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // // Compute the kernel correction fields.
  // const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  // computeCRKSPHCorrections(connectivityMap, W, mVolume, position, H, mM0, mM1, mM2, mA0, mA, mB, mC, mD, mGradA, mGradB);

  // Initialize the pressure and sound speed.
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);

  // We need to call our own evaluate derivatives method in order to initialize the 
  // viscous pressure for use in choosing a timestep.
  // vector<Physics<Dimension>*> packages(1, this);
  // State<Dimension> state(dataBase, packages);
  // StateDerivatives<Dimension> derivs(dataBase, packages);
  // this->applyGhostBoundaries(state, derivs);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();
  // this->initialize(0.0, 1.0, dataBase, state, derivs);
  // this->evaluateDerivatives(0.0, 1.0, dataBase, state, derivs);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Create the local storage for time step mask, pressure, sound speed, and correction fields.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  // dataBase.fluidPressure(mPressure);
  // dataBase.fluidSoundSpeed(mSoundSpeed);
  dataBase.resizeFluidFieldList(mPressure, 0.0,          HydroFieldNames::pressure, false);
  dataBase.resizeFluidFieldList(mSoundSpeed, 0.0,        HydroFieldNames::soundSpeed, false);
  dataBase.resizeFluidFieldList(mA,     0.0,             HydroFieldNames::A_CRKSPH, false);
  dataBase.resizeFluidFieldList(mB,     Vector::zero,    HydroFieldNames::B_CRKSPH, false);
  dataBase.resizeFluidFieldList(mGradA, Vector::zero,    HydroFieldNames::gradA_CRKSPH, false);
  dataBase.resizeFluidFieldList(mGradB, Tensor::zero,    HydroFieldNames::gradB_CRKSPH, false);
  // dataBase.resizeFluidFieldList(mSurfNorm, Vector::zero, "Surface Normal", false);

  // If we're using the compatibile energy discretization, prepare to maintain a copy
  // of the thermal energy.
  dataBase.resizeFluidFieldList(mSpecificThermalEnergy0, 0.0);
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      *mSpecificThermalEnergy0[nodeListi] = (*itr)->specificThermalEnergy();
      (*mSpecificThermalEnergy0[nodeListi]).name(HydroFieldNames::specificThermalEnergy + "0");
    }
  }

  // Now register away.
  // Mass.
  FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  state.enroll(mass);

  // Volume.
  const TableKernel<Dimension>& W = this->kernel();
  // PolicyPointer volumePolicy(new HVolumePolicy<Dimension>(W.kernelExtent()));
  // state.enroll(mVolume, volumePolicy);

  // We need to build up CompositeFieldListPolicies for the mass density and H fields
  // in order to enforce NodeList dependent limits.
  FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  boost::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > rhoPolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  boost::shared_ptr<CompositeFieldListPolicy<Dimension, SymTensor> > Hpolicy(new CompositeFieldListPolicy<Dimension, SymTensor>());
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    rhoPolicy->push_back(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                      (*itr)->rhoMax()));
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == PhysicsSpace::IntegrateH) {
      Hpolicy->push_back(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    } else {
      CHECK(HEvolution() == PhysicsSpace::IdealH);
      Hpolicy->push_back(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    }
  }
  state.enroll(massDensity, rhoPolicy);
  state.enroll(Hfield, Hpolicy);

  // Register the position update, which depends on whether we're using XSPH or not.
  FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  if (mXSPH) {
    PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
    state.enroll(position, positionPolicy);
  } else {
    PolicyPointer positionPolicy(new PositionPolicy<Dimension>());
    state.enroll(position, positionPolicy);
  }

  // Are we using the compatible energy evolution scheme?
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  if (compatibleEnergyEvolution()) {
    PolicyPointer thermalEnergyPolicy(new SpecificThermalEnergyPolicy<Dimension>(dataBase));
    // PolicyPointer thermalEnergyPolicy(new NonSymmetricSpecificThermalEnergyPolicy<Dimension>(dataBase));
    // PolicyPointer thermalEnergyPolicy(new CRKSPHSpecificThermalEnergyPolicy<Dimension>(dataBase, this->kernel()));
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,
                                                                           HydroFieldNames::specificThermalEnergy));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
    state.enroll(mSpecificThermalEnergy0);
  } else {
    PolicyPointer thermalEnergyPolicy(new IncrementFieldList<Dimension, Scalar>());
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>());
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
  }

  // Register the time step mask, initialized to 1 so that everything defaults to being
  // checked.
  state.enroll(mTimeStepMask);

  // Compute and register the pressure and sound speed.
  PolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
  PolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mSoundSpeed, csPolicy);

  // Register the CRKSPH correction fields.
  // We deliberately make these non-dynamic here.  This corrections are computed
  // during CRKSPHHydroBase::initialize, not as part of our usual state update.
  // This is necessary 'cause we need boundary conditions *and* the current set of
  // neighbors before we compute these suckers.
  state.enroll(mA);
  state.enroll(mB);
  state.enroll(mGradA);
  state.enroll(mGradB);
  // state.enroll(mSurfNorm);
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  typedef typename StateDerivatives<Dimension>::KeyType Key;
  const string DxDtName = IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const string DvDtName = IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  dataBase.resizeFluidFieldList(mEffViscousPressure, 0.0, HydroFieldNames::effectiveViscousPressure, false);
  dataBase.resizeFluidFieldList(mViscousWork, 0.0, HydroFieldNames::viscousWork, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mDmassDensityDx, Vector::zero, HydroFieldNames::massDensityGradient, false);
  dataBase.resizeFluidFieldList(mPairAccelerations, vector<Vector>(), HydroFieldNames::pairAccelerations, false);

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
  derivs.enroll(mDmassDensityDx);
  derivs.enroll(mPairAccelerations);
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

  // Compute the kernel correction fields.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  // FieldList<Dimension, Vector> surfNorm = state.fields("Surface Normal", Vector::zero);
  
  // Change CRKSPH weights here if need be!
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> vol = mass/massDensity;
  computeCRKSPHCorrections(connectivityMap, W, vol, position, H, A, B, gradA, gradB);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(A);
    (*boundItr)->applyFieldListGhostBoundary(B);
    (*boundItr)->applyFieldListGhostBoundary(gradA);
    (*boundItr)->applyFieldListGhostBoundary(gradB);
  }
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // Get the pressure and velocity gradients.
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Tensor> DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  DvDx.assignFields(CRKSPHSpace::gradientCRKSPH(velocity, position, vol, H, A, B, gradA, gradB, connectivityMap, W, NodeCoupling()));
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(DvDx);
  }

  // Get the artificial viscosity and initialize it.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();
  Q.initialize(dataBase, 
               state,
               derivs,
               this->boundaryBegin(),
               this->boundaryEnd(),
               time, 
               dt,
               W);

  // // If we're doing the RigorousSumDensity update, now is a good time to do it
  // // since we have the boundary conditions and corrections all ready to go.
  // if (densityUpdate() == PhysicsSpace::RigorousSumDensity) {
  //   const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  //   FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  //   computeCRKSPHSumMassDensity(connectivityMap, W, position, mass, vol, H, A0, A, B, massDensity);
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
  // }
  
  // // Compute the surface normal
  // surfNorm = Vector::zero;
  
  // size_t nodeListi = 0;
  // for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
  //      itr != dataBase.fluidNodeListEnd();
  //      ++itr, ++nodeListi) {
  //   const NodeList<Dimension>& nodeList = **itr;
  //   const int firstGhostNodei = nodeList.firstGhostNode();
  //   const Scalar hmin = nodeList.hmin();
  //   const Scalar hmax = nodeList.hmax();
  //   const Scalar hminratio = nodeList.hminratio();
  //   const int maxNumNeighbors = nodeList.maxNumNeighbors();
  //   const Scalar nPerh = nodeList.nodesPerSmoothingScale();
    
  //   // Iterate over the internal nodes in this NodeList.
  //   for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
  //        iItr != connectivityMap.end(nodeListi);
  //        ++iItr) {
  //     const int i = *iItr;
  //     const Vector& m1i = m1(nodeListi, i);
  //     const SymTensor& Hi = H(nodeListi, i);
  //     //const Vector& sfi = surfNorm(nodeListi, i);
      
  //     surfNorm(nodeListi,i) = Hi*m1i;
  //   }
  // }
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WQ = this->PiKernel();
  const Scalar W0 = W.kernelValue(0.0, 1.0);

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  // const FieldList<Dimension, Vector> surfNorm = state.fields("Surface Normal", Vector::zero);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(A.size() == numNodeLists);
  CHECK(B.size() == numNodeLists);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists);
  // CHECK(surfNorm.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, Vector> DrhoDx = derivatives.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, Scalar> effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  FieldList<Dimension, Scalar> viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(DrhoDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(effViscousPressure.size() == numNodeLists);
  CHECK(viscousWork.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // Size up the pair-wise accelerations before we start.
  if (mCompatibleEnergyEvolution) {
    size_t nodeListi = 0;
    for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
         itr != dataBase.fluidNodeListEnd();
         ++itr, ++nodeListi) {
      for (int i = 0; i != (*itr)->numInternalNodes(); ++i) {
        pairAccelerations(nodeListi, i).reserve(connectivityMap.numNeighborsForNode(*itr, i));
      }
    }
  }

  // Start our big loop over all FluidNodeLists.
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {
    const NodeList<Dimension>& nodeList = **itr;
    const int firstGhostNodei = nodeList.firstGhostNode();
    const Scalar hmin = nodeList.hmin();
    const Scalar hmax = nodeList.hmax();
    const Scalar hminratio = nodeList.hminratio();
    const int maxNumNeighbors = nodeList.maxNumNeighbors();
    const Scalar nPerh = nodeList.nodesPerSmoothingScale();

    // The scale for the tensile correction.
    const Scalar WnPerh = W(1.0/nPerh, 1.0);

    // Get the work field for this NodeList.
    Field<Dimension, Scalar>& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar rhoi = massDensity(nodeListi, i);
      const Scalar epsi = specificThermalEnergy(nodeListi, i);
      const Scalar Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar ci = soundSpeed(nodeListi, i);
      const Scalar Ai = A(nodeListi, i);
      const Vector& Bi = B(nodeListi, i);
      const Vector& gradAi = gradA(nodeListi, i);
      const Tensor& gradBi = gradB(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar weighti = mi/rhoi;  // Change CRKSPH weights here if need be!
      CHECK2(mi > 0.0, i << " " << mi);
      CHECK2(rhoi > 0.0, i << " " << rhoi);
      CHECK2(Ai > 0.0, i << " " << Ai);
      CHECK2(Hdeti > 0.0, i << " " << Hdeti);
      CHECK2(weighti > 0.0, i << " " << weighti);

      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      Vector& DrhoDxi = DrhoDx(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      Scalar& effViscousPressurei = effViscousPressure(nodeListi, i);
      Scalar& viscousWorki = viscousWork(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

      // // Bizarrely, in CRKSPH there is a self-contribution to gradients.  We need this 
      // // term to compute those.
      // const Scalar W0 = W.kernelValue(0.0, Hdeti);
      // const Vector selfGradContrib = W0*(Ai*Bi + gradAi);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = 1.0; // (nodeListi == nodeListj ? 1.0 : 0.2);
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();

          // Loop over the neighbors.
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;

            // Only proceed if this node pair has not been calculated yet.
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              ++ncalc;

              // Get the state for node j
              const Vector& rj = position(nodeListj, j);
              const Scalar mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar rhoj = massDensity(nodeListj, j);
              const Scalar epsj = specificThermalEnergy(nodeListj, j);
              const Scalar Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar cj = soundSpeed(nodeListj, j);
              const Scalar Aj = A(nodeListj, j);
              const Vector& Bj = B(nodeListj, j);
              const Vector& gradAj = gradA(nodeListj, j);
              const Tensor& gradBj = gradB(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar weightj = mj/rhoj;     // Change CRKSPH weights here if need be!
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Aj > 0.0 or j >= firstGhostNodej);
              CHECK(Hdetj > 0.0);
              CHECK(weightj > 0.0);

              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Vector& DrhoDxj = DrhoDx(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
              Scalar& effViscousPressurej = effViscousPressure(nodeListj, j);
              Scalar& viscousWorkj = viscousWork(nodeListj, j);
              vector<Vector>& pairAccelerationsj = pairAccelerations(nodeListj, j);
              Vector& XSPHDeltaVj = XSPHDeltaV(nodeListj, j);
              Scalar& weightedNeighborSumj = weightedNeighborSum(nodeListj, j);
              SymTensor& massSecondMomentj = massSecondMoment(nodeListj, j);

              // Node displacement.
              const Vector rij = ri - rj;
              const Vector etai = Hi*rij;
              const Vector etaj = Hj*rij;
              const Scalar etaMagi = etai.magnitude();
              const Scalar etaMagj = etaj.magnitude();
              CHECK(etaMagi >= 0.0);
              CHECK(etaMagj >= 0.0);
              const Vector vij = vi - vj;

              // Symmetrized kernel weight and gradient.
              Scalar gWi, gWj, Wi, Wj;
              Vector gradWi, gradWj;
              CRKSPHKernelAndGradient(W,  rij, -etai, Hi, Hdeti,  etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, Wj, gWj, gradWj);
              CRKSPHKernelAndGradient(W, -rij,  etaj, Hj, Hdetj, -etai, Hi, Hdeti, Aj, Bj, gradAj, gradBj, Wi, gWi, gradWi);
              const Vector deltagrad = gradWj - gradWi;
              const Vector gradWSPHi = (Hi*etai.unitVector())*WQ.gradValue(etai.magnitude(), Hdeti);
              const Vector gradWSPHj = (Hj*etaj.unitVector())*WQ.gradValue(etaj.magnitude(), Hdetj);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()/(rij2 + 1.0e-10) / FastMath::square(Dimension::pownu12(rij2 + 1.0e-10));
              weightedNeighborSumi += fweightij*std::abs(gWi);
              weightedNeighborSumj += fweightij*std::abs(gWj);
              massSecondMomenti += fweightij*gradWSPHi.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWSPHj.magnitude2()*thpt;

              // Compute the artificial viscous pressure (Pi = P/rho^2 actually).
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Scalar workQi = rhoj*rhoj*QPiij.second.dot(vij).dot(deltagrad);
              const Scalar workQj = rhoi*rhoi*QPiij.first .dot(vij).dot(deltagrad);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);
              effViscousPressurei += weightj * Qi * Wj;
              effViscousPressurej += weighti * Qj * Wi;
              viscousWorki += 0.5*weighti*weightj/mi*workQi;
              viscousWorkj += 0.5*weighti*weightj/mj*workQj;

              // Velocity gradient.
              // We've actually already set DvDx in initialize, but we need to update
              // localDvDx at this point.
              const Tensor deltaDvDxi = -weightj*vij.dyad(gradWj);
              const Tensor deltaDvDxj =  weighti*vij.dyad(gradWi);
              // DvDxi += deltaDvDxi;
              // DvDxj += deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi += deltaDvDxi;
                localDvDxj += deltaDvDxj;
              }

              // Mass density gradient.
              // DrhoDxi += weightj*rhoj*gradWj;
              // DrhoDxj += weighti*rhoi*gradWi;
              DrhoDxi += weightj*(rhoj - rhoi)*gradWj;
              DrhoDxj += weighti*(rhoi - rhoj)*gradWi;

              // // Determine an effective pressure including a term to fight the tensile instability.
              // const Scalar fij = mEpsTensile*pow(0.5*(W(etaMagi, 1.0) + W(etaMagj, 1.0))/WnPerh, mnTensile);
              // // const Scalar fij = mEpsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
              // const Scalar Ri = fij*abs(Pi);
              // const Scalar Rj = fij*abs(Pj);
              // Pi += Ri;
              // Pj += Rj;
              
              // Acceleration (CRKSPH form).
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              Vector deltaDvDti, deltaDvDtj;
              Vector forceij = 0.5*weighti*weightj*((Pi + Pj)*deltagrad + 
                                                    ((rhoi*rhoi*QPiij.first + rhoj*rhoj*QPiij.second)*deltagrad));    // <- Type III, with CRKSPH Q forces

              // // Add the filtering correction.
              // if (mfilter > 0.0) {
              //   const Vector rijhat = rij.unitVector();
              //   // const Vector nhat = sgn(rij.dot(vij))*vij.unitVector();
              //   // const Vector vijhat = vij.unitVector();
              //   // const Vector nhat = sgn(rij.dot(vij)) * vijhat;
              //   // const Scalar deltaij = min(max(0.0, max(volumeSpacing<Dimension>(mi/rhoi), volumeSpacing<Dimension>(mj/rhoj)) - rij.magnitude()),
              //   //                            vij.magnitude() * dt);
              //   const Scalar deltaij = max(0.0, 2.0*max(volumeSpacing<Dimension>(mi/rhoi), volumeSpacing<Dimension>(mj/rhoj)) - rij.magnitude());
              //   const Scalar wij = 0.5*(W.kernelValue(etaMagi, 1.0) + W.kernelValue(etaMagj, 1.0)); // /W0;
              //   const Scalar forceij0 = forceij.magnitude();
              //   forceij -= mfilter*min(forceij0, wij*mi*deltaij*2.0/(dt*dt))*rijhat;
              // }

              // // Add the filtering correction.
              // if (mfilter > 0.0) {
              //   const Scalar rijmag = rij.magnitude();
              //   const Scalar hi = rijmag/etaMagi;
              //   const Scalar hj = rijmag/etaMagj;
              //   const Scalar rknot = int(0.5*(etaMagi + etaMagj)/nPerh + 0.5) * 0.5*(hi + hj);
              //   const Scalar drij = rijmag - rknot;
              //   const Vector knothat = sgn(drij)*rij.unitVector();
              //   // const Scalar vijmag = -min(0.0, vij.dot(rijhat)); // vij.magnitude();
              //   // const Scalar ki = mi*vijmag*ci*safeInv(hi*hi);
              //   // const Scalar kj = mj*vijmag*cj*safeInv(hj*hj);
              //   const Scalar ki = mi*ci*ci*safeInv(hi*hi);
              //   const Scalar kj = mj*cj*cj*safeInv(hj*hj);
              //   // const Scalar ki = mi*vijmag*safeInv(hi*hi*ci);
              //   // const Scalar kj = mj*vijmag*safeInv(hj*hj*cj);
              //   // const Scalar drij = max(0.0, hi/nPerh - rijmag);
              //   // const Scalar drji = max(0.0, hj/nPerh - rijmag);
              //   const Scalar forceijmag = forceij.magnitude();
              //   forceij += mfilter*min(forceijmag, 0.5*(ki + kj)*abs(drij))*sgn(drij)*knothat;
              //   // if (i == 0) cerr << " --> " << j << " " << rijmag << " " << drij << " " << 0.5*mfilter*(ki + kj)*drij*knothat << endl;
              //   // forceij += min(mfilter*forceij.magnitude(), 0.5*(ki + kj))*knothat;
              // }

              deltaDvDti = -forceij/mi;
              deltaDvDtj =  forceij/mj;
              DvDti += deltaDvDti;
              DvDtj += deltaDvDtj;
              if (mCompatibleEnergyEvolution) {
                pairAccelerationsi.push_back(deltaDvDti);
                pairAccelerationsj.push_back(deltaDvDtj);
              }

              // Specific thermal energy evolution.
              DepsDti += 0.5*weighti*weightj*(Pj*vij.dot(deltagrad) + workQi)/mi;
              DepsDtj += 0.5*weighti*weightj*(Pi*vij.dot(deltagrad) + workQj)/mj;

              // Estimate of delta v (for XSPH).
              if (mXSPH and (nodeListi == nodeListj)) {
                XSPHDeltaVi -= weightj*Wj*vij;
		XSPHDeltaVj += weighti*Wi*vij;
              }

            }
          }
        }
      }
      const size_t numNeighborsi = connectivityMap.numNeighborsForNode(&nodeList, i);
      CHECK(not mCompatibleEnergyEvolution or 
            (i >= firstGhostNodei and pairAccelerationsi.size() == 0) or
            (pairAccelerationsi.size() == numNeighborsi));

      // Get the time for pairwise interactions.
      const Scalar deltaTimePair = Timing::difference(start, Timing::currentTime())/(ncalc + 1.0e-30);

      // Time evolution of the mass density.
      DrhoDti = -rhoi*DvDxi.Trace();

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        DxDti = vi + XSPHDeltaVi;
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = mSmoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             ri,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                        ri,
                                                        weightedNeighborSumi,
                                                        massSecondMomenti,
                                                        W,
                                                        hmin,
                                                        hmax,
                                                        hminratio,
                                                        nPerh,
                                                        connectivityMap,
                                                        nodeListi,
                                                        i);

      //const Scalar mag0 = DxDti.magnitude();
      const Scalar mag0 = vi.magnitude();
      //printf("MAG0=%10.3e",mag0);
      
      // if (dt > 0.0) {
      //   CHECK(m0i > 0.0);
      //   const Vector com = -m1i/m0i;
      //   const Vector dhat = com.unitVector();
      //   //const Vector delPos=com - ri;
      //   const Scalar a0 = DvDti.magnitude();
      //   const Vector delPos=com;
      //   //const Vector accel=2*delPos/(dt*dt)-2*vi/dt;
      //   //const Vector accel=2*delPos/(dt*dt)-2*DxDti/dt;
      //   const Scalar deltamag = com.magnitude();
      //   //const Vector accel=std::min(0.01*mag0, deltamag)*2*delPos/(dt*dt);
      //   const Vector accel=2*delPos/(dt*dt);
      //   const Scalar a1 = accel.magnitude();
      //   const Vector delta = mfilter*std::min(a0, a1)*accel.unitVector();
      //   //const Vector delta = 0.01*std::min(a0, a1)*accel.unitVector();

      //   // const Vector delPos2=std::min(0.01*mag0, deltamag)*dhat;
      //   //const Vector accel=2*delPos2/(dt*dt)-2*vi/dt;
      //   //const Vector accel=2*delPos/(dt*dt*mi);
      //   // printf("DVDT=%10.3e, accell=%10.3e, delta=%10.3e\n",DvDti[0],accel[0],delta[0]);
      //   // printf("COM=%10.3e, ri=%10.3e, del=%10.3e dt=%10.3e vi=%10.3e\n",com[0],ri[0],delPos[0],dt,vi[0]);
      //   //DvDti += accel;
      //   DvDti += delta;

      //   // Account for the work done as well.
      //   // DepsDti -= (vi + 0.5*dt*DvDti).dot(DvDti);
      // }

      // if (dt > 0.0) {
      //    const Vector com = centerOfMass(polyvol(nodeListi, i), DrhoDxi);
      //    const Vector dhat = com.unitVector();
      //    //const Vector delPos=com - ri;
      //    const Scalar a0 = DvDti.magnitude();
      //    const Vector delPos=com;
      //    //const Vector accel=2*delPos/(dt*dt)-2*vi/dt;
      //    //const Vector accel=2*delPos/(dt*dt)-2*DxDti/dt;
      //    const Scalar deltamag = com.magnitude();
      //    //const Vector accel=std::min(0.01*mag0, deltamag)*2*delPos/(dt*dt);
      //    const Vector accel=2*delPos/(dt*dt);
      //    const Scalar a1 = accel.magnitude();
      //    const Vector delta = mfilter*std::min(a0, a1)*accel.unitVector();
      //    //const Vector delta = 0.01*std::min(a0, a1)*accel.unitVector();

      //    // const Vector delPos2=std::min(0.01*mag0, deltamag)*dhat;
      //    //const Vector accel=2*delPos2/(dt*dt)-2*vi/dt;
      //    //const Vector accel=2*delPos/(dt*dt*mi);
      //    // printf("DVDT=%10.3e, accell=%10.3e, delta=%10.3e\n",DvDti[0],accel[0],delta[0]);
      //    // printf("COM=%10.3e, ri=%10.3e, del=%10.3e dt=%10.3e vi=%10.3e\n",com[0],ri[0],delPos[0],dt,vi[0]);
      //    //DvDti += accel;
      //    DvDti += delta;
      
      // }

      // // BLAGO!
      // if (i == 0 or i == 50) {
      //   cerr.precision(13);
      //   cerr << " --> " << i;
      //   if (i == 0) {
      //     cerr << "  ";
      //   } else {
      //     cerr << " ";
      //   }
      //   cerr << " " << Ai
      //        << " " << Bi
      //        << " " << rhoi
      //        // << " " << rhoSumi 
      //        << " " << DxDti 
      //        << " " << DrhoDti 
      //        << " " << DvDti
      //        << " " << DepsDti 
      //        << " " << DvDxi 
      //        << " " << DrhoDxi
      //        << " " << DHDti
      //        << " " << Hideali
      //        << " " << maxViscousPressurei
      //        << endl;
      // }
      // // BLAGO!

      // Increment the work for i.
      worki += Timing::difference(start, Timing::currentTime());

      // Now add the pairwise time for each neighbor we computed here.
      for (int nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const int firstGhostNodej = nodeLists[nodeListj]->firstGhostNode();
          Field<Dimension, Scalar>& workFieldj = nodeLists[nodeListj]->work();
#pragma vector always
          for (vector<int>::const_iterator jItr = connectivity.begin();
               jItr != connectivity.end();
               ++jItr) {
            const int j = *jItr;
            if (connectivityMap.calculatePairInteraction(nodeListi, i, 
                                                         nodeListj, j,
                                                         firstGhostNodej)) {
              workFieldj(j) += deltaTimePair;
            }
          }
        }
      }
    }
  }

  // BLAGO!
  // {
  //   FieldList<Dimension, Tensor> DvDx_check = gradientCRKSPH(velocity, position, vol, H, A, B, C, D, gradA, gradB, connectivityMap, W);
  //   FieldList<Dimension, Vector> gradP_check = gradientCRKSPH(pressure, position, vol, H, A, B, C, D, gradA, gradB, connectivityMap, W);
  //   for (size_t nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     for (size_t i = 0; i != DvDx[nodeListi]->numInternalElements(); ++i) {
  // 	if (std::abs(DvDx(nodeListi, i).xx() - DvDx_check(nodeListi, i).xx()) > 1e-6)
  // 	  cerr << "DVDX fail: (" << i << " " << DvDx(nodeListi, i) << " " << DvDx_check(nodeListi, i) << ") " << endl;
  // 	const Vector DvDti_check = -gradP_check(nodeListi, i)/massDensity(nodeListi, i);
  // 	if (std::abs(DvDt(nodeListi, i).x() - DvDti_check.x()) > 1.0e-6)
  // 	  cerr << "DVDT fail: (" << i << " " << DvDt(nodeListi, i) << " " << DvDti_check << ") " << endl;
  //     }
  //   }
  // }
  // BLAGO!

}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    FieldList<Dimension, Vector> accelerations = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Finalize the state after state has been updated and boundary conditions 
// enforced.  For CRKSPH this is where we update the volumes and RPKM corrections.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
postStateUpdate(const DataBase<Dimension>& dataBase,
                State<Dimension>& state,
                const StateDerivatives<Dimension>& derivs) const {

  // // Grab state we're going to use.
  // const TableKernel<Dimension>& W = this->kernel();
  // const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  // const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  // const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);

  // // Compute the volume per node.
  // FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  // computeHullVolumes(connectivityMap, position, vol);

  // // We need boundary conditions enforced on the volume before we can compute corrections.
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->applyFieldListGhostBoundary(vol);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();

  // // Compute the kernel correction fields.
  // FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  // FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  // FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CRKSPH, Vector::zero);
  // FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CRKSPH, Tensor::zero);
  // FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  // FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  // computeCRKSPHCorrections(connectivityMap, W, vol, position, H, A, B, C, D, gradA, gradB);
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) {
  //   (*boundItr)->applyFieldListGhostBoundary(A);
  //   (*boundItr)->applyFieldListGhostBoundary(B);
  //   (*boundItr)->applyFieldListGhostBoundary(C);
  //   (*boundItr)->applyFieldListGhostBoundary(D);
  //   (*boundItr)->applyFieldListGhostBoundary(gradA);
  //   (*boundItr)->applyFieldListGhostBoundary(gradB);
  // }
  // for (ConstBoundaryIterator boundItr = this->boundaryBegin();
  //      boundItr != this->boundaryEnd();
  //      ++boundItr) (*boundItr)->finalizeGhostBoundary();
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Base class finalization.
  GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (densityUpdate() == PhysicsSpace::RigorousSumDensity) {
    const TableKernel<Dimension>& W = this->kernel();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);

    // computeHullSumMassDensity(connectivityMap, this->kernel(), position, mass, H, massDensity);
    // FieldList<Dimension, Scalar> vol = mass/massDensity;
    // for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
    //      boundaryItr != this->boundaryEnd();
    //      ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(vol);
    // for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
    //      boundaryItr != this->boundaryEnd();
    //      ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
    computeCRKSPHSumMassDensity(connectivityMap, this->kernel(), position, mass, H, massDensity, NodeCoupling(), false, massDensity);

    // FieldList<Dimension, Scalar> vol = dataBase.newFluidFieldList(0.0, "volume");
    // FieldList<Dimension, FacetedVolume> polyvol = dataBase.newFluidFieldList(FacetedVolume(), "poly volume");
    // computeHullVolumes(connectivityMap, this->kernel().kernelExtent(), position, H, polyvol, vol);
    // SPHSpace::computeSPHSumMassDensity(connectivityMap, this->kernel(), position, mass, H, massDensity);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // } else if (densityUpdate() == PhysicsSpace::SumDensity) {
  //   FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   FieldList<Dimension, Scalar> massDensitySum = derivs.fields(ReplaceFieldList<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
  //                                                               HydroFieldNames::massDensity, 0.0);
  //   massDensity.assignFields(massDensitySum);
  }

  // // Add any filtering component to the node movement.
  // // Note that the FacetedVolumes are in coordinates with the node at the origin already!
  // if (mfilter > 0.0) {
  //   const TableKernel<Dimension>& W = this->kernel();
  //   const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  //   const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  //   const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  //   FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);

  //   // Find the local hulls.
  //   FieldList<Dimension, FacetedVolume> polyvol = dataBase.newFluidFieldList(FacetedVolume(), "faceted volumes");
  //   FieldList<Dimension, Scalar> vol = dataBase.newFluidFieldList(0.0, "volume");
  //   computeHullVolumes(connectivityMap, W.kernelExtent(), position, H, polyvol, vol);

  //   // Displace everyone toward their hull centroids.
  //   const unsigned numNodeLists = position.size();
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const unsigned n = position[nodeListi]->numInternalElements();
  //     for (unsigned i = 0; i != n; ++i) {
  //       const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
  //       if (mag0 > 0.0) {
  //         const Vector deltai = mfilter*polyvol(nodeListi, i).centroid();
  //         const Scalar deltamag = deltai.magnitude();
  //         const Scalar effmag = min(mag0, deltamag);
  //         position(nodeListi, i) += effmag*deltai.unitVector();
  //       }
  //     }
  //   }

  //   // Check for any boundary violations.
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
  //   this->enforceBoundaries(state, derivs);
  // }

  // // Move toward the hull neighbor centroids.
  // if (mfilter > 0.0) {
  //   const TableKernel<Dimension>& W = this->kernel();
  //   const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  //   FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  //   const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  //   const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  //   const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  //   const FieldList<Dimension, Vector> DrhoDx = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero);
  //   const unsigned numNodeLists = mass.size();
  //   const Scalar W0 = W.kernelValue(0.0, 1.0);
  //   FieldList<Dimension, Vector> delta = dataBase.newFluidFieldList(Vector::zero, "delta position");
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
  //          iItr != connectivityMap.end(nodeListi);
  //          ++iItr) {
  //       const int i = *iItr;
  //       const Vector& ri = position(nodeListi, i);
  //       const Scalar mi = mass(nodeListi, i);
  //       const Scalar rhoi = massDensity(nodeListi, i);
  //       const Vector DrhoDxi = DrhoDx(nodeListi, i);
  //       const SymTensor& Hi = H(nodeListi, i);
  //       const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
  //       const FacetedVolume polyvoli = computeNeighborHull(fullConnectivity, 1.0, ri, Hi, position);
  //       const Vector com = centerOfMass(polyvoli, DrhoDxi);
  //       delta(nodeListi, i) = mfilter*com;
  //     }
  //   }

  //   // Apply the filtering.
  //   const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const unsigned n = position[nodeListi]->numInternalElements();
  //     for (unsigned i = 0; i != n; ++i) {
  //       const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
  //       if (mag0 > 0.0) {
  //         const Scalar deltamag = delta(nodeListi, i).magnitude();
  //         const Scalar effmag = min(mfilter*mag0, deltamag);
  //         position(nodeListi, i) += effmag*delta(nodeListi, i).unitVector();
  //       }
  //     }
  //   }

  //   // Check for any boundary violations.
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
  //   this->enforceBoundaries(state, derivs);
  // }

  // This form looks for points that are too close based on specific volume.
  if (mfilter > 0.0) {
    const TableKernel<Dimension>& W = this->kernel();
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    const unsigned numNodeLists = mass.size();
    const Scalar W0 = W.kernelValue(0.0, 1.0);
    FieldList<Dimension, Vector> deltar = dataBase.newFluidFieldList(Vector::zero, "delta position");
    FieldList<Dimension, Scalar> deltav = dataBase.newFluidFieldList(0.0, "delta velocity");
    FieldList<Dimension, Scalar> weightsum = dataBase.newFluidFieldList(0.0, "delta velocity weight sum");
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
           iItr != connectivityMap.end(nodeListi);
           ++iItr) {
        const int i = *iItr;
        const Vector& ri = position(nodeListi, i);
        const Vector& vi = velocity(nodeListi, i);
        const Scalar mi = mass(nodeListi, i);
        const Scalar rhoi = massDensity(nodeListi, i);
        const SymTensor& Hi = H(nodeListi, i);
        const vector<vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
        for (unsigned nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {
          for (typename vector<int>::const_iterator jItr = fullConnectivity[nodeListj].begin();
               jItr != fullConnectivity[nodeListj].end();
               ++jItr) {
            const unsigned j = *jItr;
            const Vector& rj = position(nodeListj, j);
            const Vector& vj = velocity(nodeListj, j);
            const Scalar mj = mass(nodeListj, j);
            const Scalar rhoj = massDensity(nodeListj, j);
            const Vector rji = rj - ri;
            const Vector rjihat = rji.unitVector();
            const Scalar deltai = 2.0*max(0.0, volumeSpacing<Dimension>(mi/rhoi) + volumeSpacing<Dimension>(mj/rhoj) - rji.magnitude());
            // const Scalar deltai = max(0.0, 2.0*volumeSpacing<Dimension>((mi + mj)/(rhoi + rhoj)) - rji.magnitude());
            // deltar(nodeListi, i) -= deltai*rjihat;
            const Scalar etai = (Hi*rji).magnitude();
            const Scalar weight = W.kernelValue(etai, 1.0)/W0 * (vj - vi).magnitude();
            deltar(nodeListi, i) -= weight*deltai*rjihat;
            weightsum(nodeListi, i) += weight;
            deltav(nodeListi, i) += weight*(vj - vi).magnitude();
          }
        }
      }
    }

    // Apply the filtering.
    const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
    for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
      const unsigned n = position[nodeListi]->numInternalElements();
      for (unsigned i = 0; i != n; ++i) {
        // const Scalar hi = 1.0/(H(nodeListi, i).eigenValues().maxElement());
        // const Scalar mag0 = DvDx(nodeListi, i).eigenValues().maxAbsElement()*hi*dt;
        // const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
        const Scalar mag0 = deltav(nodeListi, i)*safeInv(weightsum(nodeListi, i))*dt;
        if (mag0 > 0.0) {
          const Scalar deltamag = deltar(nodeListi, i).magnitude();
          const Scalar effmag = mfilter*deltamag;
          // const Scalar effmag = mfilter*min(mag0, deltamag);
          position(nodeListi, i) += effmag*deltar(nodeListi, i).unitVector();
        }
      }
    }

    // Check for any boundary violations.
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
    this->enforceBoundaries(state, derivs);
  }

  // // This form looks uses Voronoi centroids.
  // if (mfilter > 0.0) {
  //   FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  //   const FieldList<Dimension, Vector> centroids = computeVoronoiCentroids(position);
  //   const FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  //   const unsigned numNodeLists = position.size();
  //   for (unsigned nodeListi = 0; nodeListi != numNodeLists; ++nodeListi) {
  //     const unsigned n = position[nodeListi]->numInternalElements();
  //     for (unsigned i = 0; i != n; ++i) {
  //       const Scalar mag0 = DxDt(nodeListi, i).magnitude() * dt;
  //       if (mag0 > 0.0) {
  //         const Vector delta = centroids(nodeListi, i) - position(nodeListi, i);
  //         const Scalar deltamag = delta.magnitude();
  //         const Scalar effmag = mfilter*min(mfilter*mag0, deltamag);
  //         position(nodeListi, i) += effmag*delta.unitVector();
  //       }
  //     }
  //   }

  //   // Check for any boundary violations.
  //   for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
  //        boundaryItr != this->boundaryEnd();
  //        ++boundaryItr) (*boundaryItr)->setAllViolationNodes(dataBase);
  //   this->enforceBoundaries(state, derivs);
  // }

}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Apply boundary conditions to the basic fluid state Fields.
  // FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);

  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  FieldList<Dimension, Vector> DrhoDx = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    // (*boundaryItr)->applyFieldListGhostBoundary(vol);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    if (compatibleEnergyEvolution()) (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
    (*boundaryItr)->applyFieldListGhostBoundary(A);
    (*boundaryItr)->applyFieldListGhostBoundary(B);
    (*boundaryItr)->applyFieldListGhostBoundary(gradA);
    (*boundaryItr)->applyFieldListGhostBoundary(gradB);
    (*boundaryItr)->applyFieldListGhostBoundary(DrhoDx);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Enforce boundary conditions on the fluid state Fields.
  // FieldList<Dimension, Scalar> vol = state.fields(HydroFieldNames::volume, 0.0);
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);

  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CRKSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CRKSPH, Vector::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CRKSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CRKSPH, Tensor::zero);
  FieldList<Dimension, Vector> DrhoDx = derivs.fields(HydroFieldNames::massDensityGradient, Vector::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    // (*boundaryItr)->enforceFieldListBoundary(vol);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    if (compatibleEnergyEvolution()) (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
    (*boundaryItr)->enforceFieldListBoundary(A);
    (*boundaryItr)->enforceFieldListBoundary(B);
    (*boundaryItr)->enforceFieldListBoundary(gradA);
    (*boundaryItr)->enforceFieldListBoundary(gradB);
    (*boundaryItr)->enforceFieldListBoundary(DrhoDx);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
dumpState(FileIO& file, string pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
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
  file.write(mDmassDensityDx, pathName + "/DmassDensityDx");
  file.write(mVolume, pathName + "/volume");
  file.write(mA, pathName + "/A");
  file.write(mB, pathName + "/B");
  file.write(mGradA, pathName + "/gradA");
  file.write(mGradB, pathName + "/gradB");
  file.write(mSurfNorm, pathName + "/surfNorm");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CRKSPHHydroBase<Dimension>::
restoreState(const FileIO& file, string pathName) {
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
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
  file.read(mDmassDensityDx, pathName + "/DmassDensityDx");
  file.read(mVolume, pathName + "/volume");
  file.read(mA, pathName + "/A");
  file.read(mB, pathName + "/B");
  file.read(mGradA, pathName + "/gradA");
  file.read(mGradB, pathName + "/gradB");
  file.read(mSurfNorm, pathName + "/surfNorm");
}

}
}

