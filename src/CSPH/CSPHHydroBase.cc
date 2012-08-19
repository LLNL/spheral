//---------------------------------Spheral++----------------------------------//
// Hydro -- The CSPH/ACSPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

#include "CSPHHydroBase.hh"
#include "CSPHUtilities.hh"
#include "computeCSPHSumMassDensity.hh"
#include "computeCSPHCorrections.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
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
#include "FileIO/FileIO.hh"

#include "TAU.h"

namespace Spheral {
namespace CSPHSpace {

using namespace std;
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

using PhysicsSpace::MassDensityType;
using PhysicsSpace::HEvolutionType;

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
CSPHHydroBase<Dimension>::
CSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
             const TableKernel<Dimension>& W,
             const TableKernel<Dimension>& WPi,
             ArtificialViscosity<Dimension>& Q,
             const double cfl,
             const bool useVelocityMagnitudeForDt,
             const bool compatibleEnergyEvolution,
             const bool XSPH,
             const MassDensityType densityUpdate,
             const HEvolutionType HUpdate):
  GenericHydro<Dimension>(W, WPi, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mSumForMassDensity(densityUpdate),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mXSPH(XSPH),
  mTimeStepMask(FieldList<Dimension, int>::Copy),
  mPressure(FieldList<Dimension, Scalar>::Copy),
  mSoundSpeed(FieldList<Dimension, Scalar>::Copy),
  mSpecificThermalEnergy0(FieldList<Dimension, Scalar>::Copy),
  mHideal(FieldList<Dimension, SymTensor>::Copy),
  mMaxViscousPressure(FieldList<Dimension, Scalar>::Copy),
  mMassDensitySum(FieldList<Dimension, Scalar>::Copy),
  mWeightedNeighborSum(FieldList<Dimension, Scalar>::Copy),
  mMassSecondMoment(FieldList<Dimension, SymTensor>::Copy),
  mXSPHDeltaV(FieldList<Dimension, Vector>::Copy),
  mDxDt(FieldList<Dimension, Vector>::Copy),
  mDvDt(FieldList<Dimension, Vector>::Copy),
  mDmassDensityDt(FieldList<Dimension, Scalar>::Copy),
  mDspecificThermalEnergyDt(FieldList<Dimension, Scalar>::Copy),
  mDHDt(FieldList<Dimension, SymTensor>::Copy),
  mDvDx(FieldList<Dimension, Tensor>::Copy),
  mInternalDvDx(FieldList<Dimension, Tensor>::Copy),
  mPairAccelerations(FieldList<Dimension, vector<Vector> >::Copy),
  mRestart(DataOutput::registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
CSPHHydroBase<Dimension>::
~CSPHHydroBase() {
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::registerState", TAU_USER);

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Create the local storage for time step mask, pressure, sound speed, and correction fields.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
  dataBase.resizeFluidFieldList(mA,     0.0,          HydroFieldNames::A_CSPH);
  dataBase.resizeFluidFieldList(mB,     Vector::zero, HydroFieldNames::B_CSPH);
  dataBase.resizeFluidFieldList(mC,     Vector::zero, HydroFieldNames::C_CSPH);
  dataBase.resizeFluidFieldList(mD,     Tensor::zero, HydroFieldNames::D_CSPH);
  dataBase.resizeFluidFieldList(mGradA, Vector::zero, HydroFieldNames::gradA_CSPH);
  dataBase.resizeFluidFieldList(mGradB, Tensor::zero, HydroFieldNames::gradB_CSPH);

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
  size_t nodeListi = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++nodeListi) {

    // Mass and mass density.
    PolicyPointer rhoPolicy(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                         (*itr)->rhoMax()));
    state.enroll((*itr)->mass());
    state.enroll((*itr)->massDensity(), rhoPolicy);

    // Register the position update, which depends on whether we're using XSPH or not.
    if (mXSPH) {
      PolicyPointer positionPolicy(new IncrementState<Dimension, Vector>());
      state.enroll((*itr)->positions(), positionPolicy);
    } else {
      PolicyPointer positionPolicy(new PositionPolicy<Dimension>());
      state.enroll((*itr)->positions(), positionPolicy);
    }

    // Are we using the compatible energy evolution scheme?
    if (compatibleEnergyEvolution()) {
      PolicyPointer thermalEnergyPolicy(new SpecificThermalEnergyPolicy<Dimension>(dataBase));
      PolicyPointer velocityPolicy(new IncrementState<Dimension, Vector>(HydroFieldNames::position,
                                                                         HydroFieldNames::specificThermalEnergy));

      state.enroll((*itr)->specificThermalEnergy(), thermalEnergyPolicy);
      state.enroll((*itr)->velocity(), velocityPolicy);
      state.enroll(*mSpecificThermalEnergy0[nodeListi]);
    } else {
      PolicyPointer thermalEnergyPolicy(new IncrementState<Dimension, Scalar>());
      PolicyPointer velocityPolicy(new IncrementState<Dimension, Vector>());
      state.enroll((*itr)->specificThermalEnergy(), thermalEnergyPolicy);
      state.enroll((*itr)->velocity(), velocityPolicy);
    }

    // Register the H tensor.
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == PhysicsSpace::IntegrateH) {
      PolicyPointer Hpolicy(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
      state.enroll((*itr)->Hfield(), Hpolicy);
    } else {
      CHECK(HEvolution() == PhysicsSpace::IdealH);
      PolicyPointer Hpolicy(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
      state.enroll((*itr)->Hfield(), Hpolicy);
    }

    // Register the time step mask, initialized to 1 so that everything defaults to being
    // checked.
    state.enroll(*mTimeStepMask[nodeListi]);

    // Compute and register the pressure and sound speed.
    PolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
    PolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
    state.enroll(*mPressure[nodeListi], pressurePolicy);
    state.enroll(*mSoundSpeed[nodeListi], csPolicy);

    // Register the CSPH correction fields.
    // We deliberately make these non-dynamic here.  This corrections are computed
    // during CSPHHydroBase::initialize, not as part of our usual state update.
    state.enroll(*mA[nodeListi]);
    state.enroll(*mB[nodeListi]);
    state.enroll(*mC[nodeListi]);
    state.enroll(*mD[nodeListi]);
    state.enroll(*mGradA[nodeListi]);
    state.enroll(*mGradB[nodeListi]);
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::registerDerivatives", TAU_USER);

  typedef typename StateDerivatives<Dimension>::KeyType Key;
  const string DxDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::position;
  const string DvDtName = IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity;

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  // the ArtificialVisocisity::initialize step).
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  dataBase.resizeFluidFieldList(mMassDensitySum, 0.0, ReplaceState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mPairAccelerations, vector<Vector>(), HydroFieldNames::pairAccelerations, false);

  size_t i = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++i) {
    derivs.enroll(*mHideal[i]);
    derivs.enroll(*mMaxViscousPressure[i]);
    derivs.enroll(*mMassDensitySum[i]);
    derivs.enroll(*mWeightedNeighborSum[i]);
    derivs.enroll(*mMassSecondMoment[i]);
    derivs.enroll(*mXSPHDeltaV[i]);

    // These two (the position and velocity updates) may be registered
    // by other physics packages as well, so we need to be careful
    // not to duplicate if so.
    const Key DxDtKey = State<Dimension>::buildFieldKey(DxDtName, (*itr)->name());
    const Key DvDtKey = State<Dimension>::buildFieldKey(DvDtName, (*itr)->name());
    if (not derivs.registered(DxDtKey)) derivs.enroll(*mDxDt[i]);
    if (not derivs.registered(DvDtKey)) derivs.enroll(*mDvDt[i]);

    derivs.enroll(*mDmassDensityDt[i]);
    derivs.enroll(*mDspecificThermalEnergyDt[i]);
    derivs.enroll(*mDHDt[i]);
    derivs.enroll(*mDvDx[i]);
    derivs.enroll(*mInternalDvDx[i]);
    derivs.enroll(*mPairAccelerations[i]);
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::initialize", TAU_USER);

  // Compute the kernel correction fields.
  const TableKernel<Dimension>& W = this->kernel();
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);
  computeCSPHCorrections(connectivityMap, W, mass, position, H, A, B, C, D, gradA, gradB);
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr) {
    (*boundItr)->applyFieldListGhostBoundary(A);
    (*boundItr)->applyFieldListGhostBoundary(B);
    (*boundItr)->applyFieldListGhostBoundary(C);
    (*boundItr)->applyFieldListGhostBoundary(D);
    (*boundItr)->applyFieldListGhostBoundary(gradA);
    (*boundItr)->applyFieldListGhostBoundary(gradB);
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

  // We depend on the caller knowing to finalize the ghost boundaries!
}

//------------------------------------------------------------------------------
// Determine the principle derivatives for the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::evaluateDerivatives", TAU_USER);
  TAU_PROFILE_TIMER(TimeCSPHHydroBaseCalcDerivs, "CSPHHydroBase", "::evaluateDerivatives : calc derivatives", TAU_USER);
  TAU_PROFILE_TIMER(TimeCSPHHydroBaseGetFieldLists, "CSPHNode", "::evaluateDerivatives : get FieldLists", TAU_USER);
  TAU_PROFILE_TIMER(TimeCSPHHydroBaseNodeIState, "CSPHNode", "::evaluateDerivatives : node i state", TAU_USER);

  // Get the ArtificialViscosity.
  ArtificialViscosity<Dimension>& Q = this->artificialViscosity();

  // The kernels and such.
  const TableKernel<Dimension>& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  TAU_PROFILE_START(TimeCSPHHydroBaseGetFieldLists);
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
  const FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  const FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  const FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  const FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  const FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  const FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);
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
  CHECK(C.size() == numNodeLists);
  CHECK(D.size() == numNodeLists);
  CHECK(gradA.size() == numNodeLists);
  CHECK(gradB.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.fields(ReplaceState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DxDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Vector> DvDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> DepsDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Tensor> DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  FieldList<Dimension, Tensor> localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  FieldList<Dimension, SymTensor> DHDt = derivatives.fields(IncrementState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, SymTensor> Hideal = derivatives.fields(ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, SymTensor::zero);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.fields(HydroFieldNames::pairAccelerations, vector<Vector>());
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  CHECK(rhoSum.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(maxViscousPressure.size() == numNodeLists);
  CHECK(pairAccelerations.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  TAU_PROFILE_STOP(TimeCSPHHydroBaseGetFieldLists);

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
  TAU_PROFILE_START(TimeCSPHHydroBaseCalcDerivs);
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

    // Get the work field for this NodeList.
    Field<Dimension, Scalar>& workFieldi = nodeList.work();

    // Iterate over the internal nodes in this NodeList.
    for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
         iItr != connectivityMap.end(nodeListi);
         ++iItr) {
      const int i = *iItr;

      TAU_PROFILE_START(TimeCSPHHydroBaseNodeIState);
      // Prepare to accumulate the time.
      const Time start = Timing::currentTime();
      size_t ncalc = 0;

      // Get the state for node i.
      const Vector& ri = position(nodeListi, i);
      const Scalar& mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Scalar& epsi = specificThermalEnergy(nodeListi, i);
      const Scalar& Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar& ci = soundSpeed(nodeListi, i);
      const Scalar& Ai = A(nodeListi, i);
      const Vector& Bi = B(nodeListi, i);
      const Vector& gradAi = gradA(nodeListi, i);
      const Tensor& gradBi = gradB(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar weighti = mi;
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Ai > 0.0);
      CHECK(Hdeti > 0.0);

      Scalar& rhoSumi = rhoSum(nodeListi, i);
      Vector& DxDti = DxDt(nodeListi, i);
      Scalar& DrhoDti = DrhoDt(nodeListi, i);
      Vector& DvDti = DvDt(nodeListi, i);
      Scalar& DepsDti = DepsDt(nodeListi, i);
      Tensor& DvDxi = DvDx(nodeListi, i);
      Tensor& localDvDxi = localDvDx(nodeListi, i);
      SymTensor& DHDti = DHDt(nodeListi, i);
      SymTensor& Hideali = Hideal(nodeListi, i);
      Scalar& maxViscousPressurei = maxViscousPressure(nodeListi, i);
      vector<Vector>& pairAccelerationsi = pairAccelerations(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);
      TAU_PROFILE_STOP(TimeCSPHHydroBaseNodeIState);

      // Iterate over the NodeLists.
      for (size_t nodeListj = 0; nodeListj != numNodeLists; ++nodeListj) {

        // Connectivity of this node with this NodeList.  We only need to proceed if
        // there are some nodes in this list.
        const vector<int>& connectivity = fullConnectivity[nodeListj];
        if (connectivity.size() > 0) {
          const double fweightij = (nodeListi == nodeListj ? 1.0 : 0.2);
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
              const Scalar& mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar& rhoj = massDensity(nodeListj, j);
              const Scalar& epsj = specificThermalEnergy(nodeListj, j);
              const Scalar& Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar& cj = soundSpeed(nodeListj, j);
              const Scalar& Aj = A(nodeListj, j);
              const Vector& Bj = B(nodeListj, j);
              const Vector& gradAj = gradA(nodeListj, j);
              const Tensor& gradBj = gradB(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar weightj = mj;
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
              CHECK(Aj > 0.0 or j >= firstGhostNodej);
              CHECK(Hdetj > 0.0);

              Scalar& rhoSumj = rhoSum(nodeListj, j);
              Vector& DxDtj = DxDt(nodeListj, j);
              Scalar& DrhoDtj = DrhoDt(nodeListj, j);
              Vector& DvDtj = DvDt(nodeListj, j);
              Scalar& DepsDtj = DepsDt(nodeListj, j);
              Tensor& DvDxj = DvDx(nodeListj, j);
              Tensor& localDvDxj = localDvDx(nodeListj, j);
              Scalar& maxViscousPressurej = maxViscousPressure(nodeListj, j);
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

              // Symmetrized kernel weight and gradient.
              Scalar Wi, gWi, Wj, gWj;
              Vector gradWi, gradWj;
              CSPHKernelAndGradient(W, rij, etaj, Hj, Hdetj, Ai, Bi, gradAi, gradBi, Wj, gWj, gradWj);
              CSPHKernelAndGradient(W, rij, etai, Hi, Hdeti, Aj, Bj, gradAj, gradBj, Wi, gWi, gradWi);

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()/(rij2 + 1.0e-10) / FastMath::square(Dimension::pownu12(rij2 + 1.0e-10));
              const Vector gradWiSPH = Hi*etai.unitVector() * gWi;
              const Vector gradWjSPH = Hj*etaj.unitVector() * gWj;
              weightedNeighborSumi += fweightij*std::abs(gWi);
              weightedNeighborSumj += fweightij*std::abs(gWj);
              massSecondMomenti += fweightij*gradWiSPH.magnitude2()*thpt;
              massSecondMomentj += fweightij*gradWjSPH.magnitude2()*thpt;

              // Contribution to the sum density (only if the same material).
              if (nodeListi == nodeListj) {
                rhoSumi += mj*W.kernelValue(etaMagi, Hdeti);
                rhoSumj += mi*W.kernelValue(etaMagj, Hdetj);
              }

              // Mass density evolution.
              const Vector vij = vi - vj;
              const double deltaDrhoDti = weightj*rhoj*vij.dot(gradWj);
              const double deltaDrhoDtj = weighti*rhoi*vij.dot(gradWi);
              DrhoDti += deltaDrhoDti;
              DrhoDtj += deltaDrhoDtj;

              // Compute the pair-wise artificial viscosity.
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Vector Qacci = 0.5*(QPiij.first *gradWi);
              const Vector Qaccj = 0.5*(QPiij.second*gradWj);
              const Scalar workQi = 0.5*(QPiij.first *vij).dot(gradWi);
              const Scalar workQj = 0.5*(QPiij.second*vij).dot(gradWj);
              // const Scalar workQi = vij.dot(Qacci);
              // const Scalar workQj = vij.dot(Qaccj);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);

              // Acceleration.
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const Vector deltaDvDti = weightj*rhoj*(Pi - Pj)*gradWj/(rhoi*rhoi) + Qacci + Qaccj;
              const Vector deltaDvDtj = weighti*rhoi*(Pi - Pj)*gradWi/(rhoj*rhoj) + Qacci + Qaccj;
              DvDti += deltaDvDti;
              DvDtj += deltaDvDtj;

              // Specific thermal energy evolution.
              if (mCompatibleEnergyEvolution) {
                if (i < firstGhostNodei) pairAccelerationsi.push_back(deltaDvDti);
                if (j < firstGhostNodej) pairAccelerationsj.push_back(deltaDvDtj);
              }

              // Velocity gradient.
              const Tensor deltaDvDxi = weightj*rhoj*vij.dyad(gradWj);
              const Tensor deltaDvDxj = weighti*rhoi*vij.dyad(gradWi);
              DvDxi -= deltaDvDxi;
              DvDxj -= deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= deltaDvDxi;
                localDvDxj -= deltaDvDxj;
              }

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

      // Finish the mass density sum.
      rhoSumi = Ai*(rhoSumi + mi*W(0.0, Hdeti));

      // Finish the velocity gradient.
      CHECK(rhoi > 0.0);
      DvDxi /= rhoi;
      localDvDxi /= rhoi;

      // The specific thermal energy evolution.
//     DepsDti = Pi/(rhoi*rhoi)*DrhoDti;
      DepsDti = -Pi/rhoi * DvDxi.Trace();

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
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh,
                                                             maxNumNeighbors);
      Hideali = mSmoothingScaleMethod.newSmoothingScale(Hi,
                                                        weightedNeighborSumi,
                                                        massSecondMomenti,
                                                        numNeighborsi,
                                                        W,
                                                        hmin,
                                                        hmax,
                                                        hminratio,
                                                        nPerh,
                                                        maxNumNeighbors);

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
  TAU_PROFILE_STOP(TimeCSPHHydroBaseCalcDerivs);
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::finalizeDerivatives", TAU_USER);

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    FieldList<Dimension, Vector> accelerations = derivs.fields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
}

//------------------------------------------------------------------------------
// Finalize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::finalize", TAU_USER);

  // Base class finalization.
  GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (sumForMassDensity() == PhysicsSpace::RigorousSumDensity) {
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Vector> position = state.fields(HydroFieldNames::position, Vector::zero);
    const FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
    const FieldList<Dimension, SymTensor> H = state.fields(HydroFieldNames::H, SymTensor::zero);
    const FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    computeCSPHSumMassDensity(connectivityMap, this->kernel(), position, mass, H, A, massDensity);
  } else if (sumForMassDensity() == PhysicsSpace::SumDensity) {
    FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    FieldList<Dimension, Scalar> massDensitySum = derivs.fields(ReplaceState<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                                HydroFieldNames::massDensity, 0.0);
    massDensity.assignFields(massDensitySum);
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::applyGhostBoundaries", TAU_USER);

  // Apply boundary conditions to the basic fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);

  FieldList<Dimension, Scalar> A = state.fields(HydroFieldNames::A_CSPH, 0.0);
  FieldList<Dimension, Vector> B = state.fields(HydroFieldNames::B_CSPH, Vector::zero);
  FieldList<Dimension, Vector> C = state.fields(HydroFieldNames::C_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> D = state.fields(HydroFieldNames::D_CSPH, Tensor::zero);
  FieldList<Dimension, Vector> gradA = state.fields(HydroFieldNames::gradA_CSPH, Vector::zero);
  FieldList<Dimension, Tensor> gradB = state.fields(HydroFieldNames::gradB_CSPH, Tensor::zero);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    if (compatibleEnergyEvolution()) (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
    (*boundaryItr)->applyFieldListGhostBoundary(A);
    (*boundaryItr)->applyFieldListGhostBoundary(B);
    (*boundaryItr)->applyFieldListGhostBoundary(C);
    (*boundaryItr)->applyFieldListGhostBoundary(D);
    (*boundaryItr)->applyFieldListGhostBoundary(gradA);
    (*boundaryItr)->applyFieldListGhostBoundary(gradB);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // TAU timers.
  TAU_PROFILE("CSPHHydroBase", "::enforceBoundaries", TAU_USER);

  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.fields(HydroFieldNames::specificThermalEnergy + "0", 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    if (compatibleEnergyEvolution()) (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CSPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
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
CSPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
}

}
}

