//---------------------------------Spheral++----------------------------------//
// TotalHydro -- The Hydro class for methods that evolve total conserved
// quantities (mass, momentum, and energy) rather than specific quantities.
//
// Created by JMO, Mon Jul 12 21:07:52 PDT 2010
//----------------------------------------------------------------------------//
#include "TotalHydro.hh"
#include "HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/NonDynamicState.hh"
#include "SpecificThermalEnergyPolicy.hh"
#include "PositionPolicy.hh"
#include "PressurePolicy.hh"
#include "SoundSpeedPolicy.hh"
#include "ArtificialViscosity/ArtificialViscosity.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/iterateIdealH.hh"
#include "FileIO/FileIO.hh"

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>

namespace Spheral {


//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
TotalHydro<Dimension>::
TotalHydro(const TableKernel<Dimension>& W,
           const TableKernel<Dimension>& WPi,
           ArtificialViscosity<Dimension>& Q,
           const HEvolutionType HUpdate,
           const double hmin,
           const double hmax,
           const double hratiomin):
  GenericHydro<Dimension>(W, WPi, Q),
  mHEvolution(HUpdate),
  mhmin(hmin),
  mhmax(hmax),
  mhratiomin(hratiomin),
  mHideal(FieldList<Dimension, SymTensor>::Copy),
  mTimeStepMask(FieldList<Dimension, int>::Copy),
  mPressure(FieldList<Dimension, Scalar>::Copy),
  mSoundSpeed(FieldList<Dimension, Scalar>::Copy),
  mPositionWeight(FieldList<Dimension, Scalar>::Copy),
  mWeightedNeighborSum(FieldList<Dimension, Scalar>::Copy),
  mVolume(FieldList<Dimension, Scalar>::Copy),
  mTotalEnergy(FieldList<Dimension, Scalar>::Copy),
  mDVDt(FieldList<Dimension, Scalar>::Copy),
  mDEDt(FieldList<Dimension, Scalar>::Copy),
  mLinearMomentum(FieldList<Dimension, Vector>::Copy),
  mDpmomDt(FieldList<Dimension, Vector>::Copy),
  mMassSecondMoment(FieldList<Dimension, SymTensor>::Copy),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
TotalHydro<Dimension>::
~TotalHydro() {
}

//------------------------------------------------------------------------------
// Determine the principle derivatives for the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {

  REQUIRE(this->valid());

  // Iterate over the FluidNodeLists.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {

    // Have the FluidNodeList set it's own derivatives.
    (*itr)->calculateDerivatives(time, 
                                 dt,
                                 connectivityMap,
                                 state,
                                 derivatives);

  }
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::IntPolicyPointerType IntPolicyPointer;
  typedef typename State<Dimension>::ScalarPolicyPointerType ScalarPolicyPointer;
  typedef typename State<Dimension>::VectorPolicyPointerType VectorPolicyPointer;
  typedef typename State<Dimension>::TensorPolicyPointerType TensorPolicyPointer;
  typedef typename State<Dimension>::SymTensorPolicyPointerType SymTensorPolicyPointer;

  // Compute the total momentum and energy
  dataBase.resizeFluidFieldList(mVolume, 0.0, HydroFieldNames::totalEnergy);
  dataBase.resizeFluidFieldList(mLinearMomentum, Vector::zero, HydroFieldNames::linearMomentum);
  dataBase.resizeFluidFieldList(mTotalEnergy, 0.0, HydroFieldNames::totalEnergy);

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  mPressure = FieldList<Dimension, Scalar>(FieldList<Dimension, Scalar>::Copy);
  mSoundSpeed = FieldList<Dimension, Scalar>(FieldList<Dimension, Scalar>::Copy);
  dataBase.resizeFluidFieldList(mPositionWeight, 1.0, HydroFieldNames::positionWeight);

  // Build the min and max H tensors.
  const double hmin = this->hmin();
  const double hmax = this->hmax();
  CHECK(distinctlyGreaterThan(hmin, 0.0) &&
        distinctlyGreaterThan(hmax, 0.0) &&
        hmax >= hmin);
  const Scalar hminInv = 1.0/hmin;
  const Scalar hmaxInv = 1.0/hmax;
  CHECK(hminInv >= hmaxInv);

  // Now register away.
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    Field<Dimension, Scalar>& m = (*itr)->mass();
    const Field<Dimension, Scalar>& rho = (*itr)->massDensity();
    const Field<Dimension, Vector>& vel = (*itr)->velocity();
    const Field<Dimension, Scalar>& eps = (*itr)->specificThermalEnergy();

    // Compute the volume, momentum, and energy.
    typename FieldList<Dimension, Scalar>::iterator VItr = mVolume.fieldForNodeList(**itr);
    typename FieldList<Dimension, Vector>::iterator pmomItr = mLinearMomentum.fieldForNodeList(**itr);
    typename FieldList<Dimension, Scalar>::iterator EItr = mTotalEnergy.fieldForNodeList(**itr);
    CHECK(VItr != mVolume.end());
    CHECK(pmomItr != mLinearMomentum.end());
    CHECK(EItr != mTotalEnergy.end());
    Field<Dimension, Scalar>& V = **VItr;
    Field<Dimension, Vector>& pmom = **pmomItr;
    Field<Dimension, Scalar>& E = **EItr;
    for (size_t i = 0; i != (*itr)->numNodes(); ++i) {
      V[i] = m[i]/rho[i];
      pmom[i] = m[i]*vel[i];
      E[i] = m[i]*(0.5*vel[i].magnitude2() + eps[i]);
    }

    // Mass.
    ScalarPolicyPointer massPolicy(new NonDynamicState<Dimension, Scalar>());
    state.registerField(m, massPolicy);

    // Register the position update.
    VectorPolicyPointer positionPolicy(new IncrementState<Dimension, Vector>());
    state.registerField((*itr)->positions(), positionPolicy);

    // Register the volume.
    ScalarPolicyPointer VPolicy(new IncrementState<Dimension, Scalar>());
    state.registerField(V, VPolicy);

    // Register the momentum.
    VectorPolicyPointer pmomPolicy(new IncrementState<Dimension, Vector>());
    state.registerField(pmom, pmomPolicy);

    // Register the total energy.
    ScalarPolicyPointer EPolicy(new IncrementState<Dimension, Scalar>());
    state.registerField(E, EPolicy);

    // Register the H tensor.
    if (HEvolution() == IntegrateH) {
      SymTensorPolicyPointer Hpolicy(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
      state.registerField((*itr)->Hfield(), Hpolicy);
    } else {
      CHECK(HEvolution() == IdealH);
      SymTensorPolicyPointer Hpolicy(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
      state.registerField((*itr)->Hfield(), Hpolicy);
    }

    // Register the time step mask, initialized to 1 so that everything defaults to being
    // checked.
    IntPolicyPointer timeStepMaskPolicy(new NonDynamicState<Dimension, int>());
    typename FieldList<Dimension, int>::iterator maskItr = mTimeStepMask.fieldForNodeList(**itr);
    CHECK(maskItr != mTimeStepMask.end());
    state.registerField(**maskItr, timeStepMaskPolicy);

    // Compute and register the pressure and sound speed.
    ScalarPolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
    ScalarPolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
    mPressure.appendField((*itr)->pressure());
    mSoundSpeed.appendField((*itr)->soundSpeed());
    typename FieldList<Dimension, Scalar>::iterator PItr = mPressure.fieldForNodeList(**itr);
    typename FieldList<Dimension, Scalar>::iterator csItr = mSoundSpeed.fieldForNodeList(**itr);
    CHECK(PItr < mPressure.end());
    CHECK(csItr < mSoundSpeed.end());
    state.registerField(**PItr, pressurePolicy);
    state.registerField(**csItr, csPolicy);

    // Register the position weight for the H measurement.
    ScalarPolicyPointer positionWeightPolicy(new NonDynamicState<Dimension, Scalar>());
    typename FieldList<Dimension, Scalar>::iterator positionWeightItr = mPositionWeight.fieldForNodeList(**itr);
    CHECK(positionWeightItr < mPositionWeight.end());
    state.registerField(**positionWeightItr, positionWeightPolicy);
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  // Create the scratch fields.
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment);
  dataBase.resizeFluidFieldList(mDVDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::volume);
  dataBase.resizeFluidFieldList(mDpmomDt, Vector::zero, IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::linearMomentum);
  dataBase.resizeFluidFieldList(mDEDt, 0.0, IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::totalEnergy);

  size_t i = 0;
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr, ++i) {
    derivs.registerField((*itr)->DxDt());
    derivs.registerField((*itr)->DvelocityDx());
    derivs.registerField((*itr)->internalDvDx());
    derivs.registerField((*itr)->DHDt());
    derivs.registerField(*mDVDt[i]);
    derivs.registerField(*mDpmomDt[i]);
    derivs.registerField(*mDEDt[i]);
    derivs.registerField(*mHideal[i]);
    derivs.registerField(*mWeightedNeighborSum[i]);
    derivs.registerField(*mMassSecondMoment[i]);
  }

}

//------------------------------------------------------------------------------
// Post-state update.  Fill in the NodeList specific state based on the total
// conserved quantities.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
postStateUpdate(const Scalar time, 
                const Scalar dt,
                const DataBase<Dimension>& dataBase, 
                State<Dimension>& state,
                StateDerivatives<Dimension>& derivatives) {

  // Walk the FluidNodeLists.
  for (typename DataBase<Dimension>::ConstFluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    const Field<Dimension, Scalar>& m = (*itr)->mass();
    Field<Dimension, Scalar>& rho = const_cast<Field<Dimension, Scalar>&>((*itr)->massDensity());
    Field<Dimension, Vector>& vel = const_cast<Field<Dimension, Vector>&>((*itr)->velocity());
    Field<Dimension, Scalar>& eps = const_cast<Field<Dimension, Scalar>&>((*itr)->specificThermalEnergy());

    // Extract the updated momentum and energy from the state.
    typename FieldList<Dimension, Scalar>::iterator VItr = mVolume.fieldForNodeList(**itr);
    typename FieldList<Dimension, Vector>::iterator pmomItr = mLinearMomentum.fieldForNodeList(**itr);
    typename FieldList<Dimension, Scalar>::iterator EItr = mTotalEnergy.fieldForNodeList(**itr);
    CHECK(VItr != mVolume.end());
    CHECK(pmomItr != mLinearMomentum.end());
    CHECK(EItr != mTotalEnergy.end());
    const Field<Dimension, Scalar>& V = **VItr;
    const Field<Dimension, Vector>& pmom = **pmomItr;
    const Field<Dimension, Scalar>& E = **EItr;

    // Now update the velocity and specific thermal energy.
    for (size_t i = 0; i != (*itr)->numInternalNodes(); ++i) {
      rho[i] = m[i]/V[i];
      vel[i] = pmom[i]/m[i];
      eps[i] = E[i]/m[i] - 0.5*vel[i].magnitude2();
    }
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  REQUIRE(this->valid());

  // Apply boundary conditions to the basic fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  FieldList<Dimension, Scalar> V = state.scalarFields(HydroFieldNames::volume);
  FieldList<Dimension, Vector> pmom = state.vectorFields(HydroFieldNames::linearMomentum);
  FieldList<Dimension, Scalar> E = state.scalarFields(HydroFieldNames::totalEnergy);
  FieldList<Dimension, Scalar> pressure = state.scalarFields(HydroFieldNames::pressure);
  FieldList<Dimension, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);
  FieldList<Dimension, Scalar> positionWeight = state.scalarFields(HydroFieldNames::positionWeight);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(V);
    (*boundaryItr)->applyFieldListGhostBoundary(pmom);
    (*boundaryItr)->applyFieldListGhostBoundary(E);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(positionWeight);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  REQUIRE(this->valid());

  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  FieldList<Dimension, Scalar> V = state.scalarFields(HydroFieldNames::volume);
  FieldList<Dimension, Vector> pmom = state.vectorFields(HydroFieldNames::linearMomentum);
  FieldList<Dimension, Scalar> E = state.scalarFields(HydroFieldNames::totalEnergy);
  FieldList<Dimension, Scalar> pressure = state.scalarFields(HydroFieldNames::pressure);
  FieldList<Dimension, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(V);
    (*boundaryItr)->enforceFieldListBoundary(pmom);
    (*boundaryItr)->enforceFieldListBoundary(E);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
}  

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
TotalHydro<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
}

}

//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class TotalHydro< Dim<1> >;
  template class TotalHydro< Dim<2> >;
  template class TotalHydro< Dim<3> >;
}
