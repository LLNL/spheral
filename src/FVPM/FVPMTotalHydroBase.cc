//---------------------------------Spheral++----------------------------------//
// FVPMTotalHydroBase -- The FVPM/AFVPM hydrodynamic package for Spheral++.
//
// Created by JMO, Tue Aug  3 21:22:10 PDT 2010
//----------------------------------------------------------------------------//
#include "FVPMTotalHydroBase.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/GenericHydro.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/NonDynamicState.hh"
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
FVPMTotalHydroBase<Dimension>::
FVPMTotalHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
                   const TableKernel<Dimension>& W,
                   ArtificialViscosity<Dimension>& Q,
                   const double cfl,
                   const bool useVelocityMagnitudeForDt,
                   const bool compatibleEnergyEvolution,
                   const MassDensityType densityUpdate):
  GenericHydro<Dimension>(W, W, Q, cfl, useVelocityMagnitudeForDt),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mTimeStepMask(FieldList<Dimension, int>::Copy),
  mPressure(FieldList<Dimension, Scalar>::Copy),
  mSoundSpeed(FieldList<Dimension, Scalar>::Copy),
  mPositionWeight(FieldList<Dimension, Scalar>::Copy),
  mOmegaGradh(FieldList<Dimension, Scalar>::Copy),
  mSpecificThermalEnergy0(FieldList<Dimension, Scalar>::Copy),
  mHideal(FieldList<Dimension, SymTensor>::Copy),
  mMaxViscousPressure(FieldList<Dimension, Scalar>::Copy),
  mMassDensitySum(FieldList<Dimension, Scalar>::Copy),
  mWeightedNeighborSum(FieldList<Dimension, Scalar>::Copy),
  mMassSecondMoment(FieldList<Dimension, SymTensor>::Copy),
  mXSPHWeightSum(FieldList<Dimension, Scalar>::Copy),
  mXSPHDeltaV(FieldList<Dimension, Vector>::Copy),
  mDxDt(FieldList<Dimension, Vector>::Copy),
  mDvDt(FieldList<Dimension, Vector>::Copy),
  mDmassDensityDt(FieldList<Dimension, Scalar>::Copy),
  mDspecificThermalEnergyDt(FieldList<Dimension, Scalar>::Copy),
  mDHDt(FieldList<Dimension, SymTensor>::Copy),
  mDvDx(FieldList<Dimension, Tensor>::Copy),
  mInternalDvDx(FieldList<Dimension, Tensor>::Copy),
  mPairAccelerations(FieldList<Dimension, vector<Vector> >::Copy),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
FVPMTotalHydroBase<Dimension>::
~FVPMTotalHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // Create storage for the pressure and sound speed.
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);

  // Initialize the pressure and sound speed.
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
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
  dataBase.fluidVolume(mVolume);
  dataBase.fluidLinearMomentum(mLinearMomentum);
  dataBase.fluidEnergy(mTotalEnergy);

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
  dataBase.resizeFluidFieldList(mPositionWeight, 1.0, HydroFieldNames::positionWeight);

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

    // Mass and volume.
    ScalarPolicyPointer massPolicy(new NonDynamicState<Dimension, Scalar>());
    ScalarPolicyPointer VPolicy(new IncrementState<Dimension, Scalar>());
    state.registerField((*itr)->mass(), massPolicy);
    state.registerField(*mVolume[nodeListi], VPolicy);

    // Position and momentum.
    VectorPolicyPointer positionPolicy(new IncrementState<Dimension, Vector>());
    VectorPolicyPointer pmomPolicy(new IncrementState<Dimension, Vector>());
    state.registerField((*itr)->positions(), positionPolicy);
    state.registerField(*mLinearMomentum[nodeListi], pmomPolicy);

    // Register the total energy.
    ScalarPolicyPointer EPolicy(new IncrementState<Dimension, Scalar>());
    state.registerField(*mTotalEnergy[nodeListi], EPolicy);

    // Register the H tensor.
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
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
    state.registerField(*mTimeStepMask[nodeListi], timeStepMaskPolicy);

    // Compute and register the pressure and sound speed.
    ScalarPolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
    ScalarPolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
    state.registerField(*mPressure[nodeListi], pressurePolicy);
    state.registerField(*mSoundSpeed[nodeListi], csPolicy);

    // Register the position weight for the H measurement.
    ScalarPolicyPointer positionWeightPolicy(new NonDynamicState<Dimension, Scalar>());
    state.registerField(*mPositionWeight[nodeListi], positionWeightPolicy);
  }
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {

  typedef typename StateDerivatives<Dimension>::FieldKeyType Key;
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
  dataBase.resizeFluidFieldList(mXSPHWeightSum, 0.0, HydroFieldNames::XSPHWeightSum, false);
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
    derivs.registerField(*mHideal[i]);
    derivs.registerField(*mMaxViscousPressure[i]);
    derivs.registerField(*mMassDensitySum[i]);
    derivs.registerField(*mWeightedNeighborSum[i]);
    derivs.registerField(*mMassSecondMoment[i]);
    derivs.registerField(*mXSPHWeightSum[i]);
    derivs.registerField(*mXSPHDeltaV[i]);

    // These two (the position and velocity updates) may be registered
    // by other physics packages as well, so we need to be careful
    // not to duplicate if so.
    const Key DxDtKey(*itr, DxDtName);
    const Key DvDtKey(*itr, DvDtName);
    if (not derivs.vectorFieldRegistered(DxDtKey)) derivs.registerField(*mDxDt[i]);
    if (not derivs.vectorFieldRegistered(DvDtKey)) derivs.registerField(*mDvDt[i]);

    derivs.registerField(*mDmassDensityDt[i]);
    derivs.registerField(*mDspecificThermalEnergyDt[i]);
    derivs.registerField(*mDHDt[i]);
    derivs.registerField(*mDvDx[i]);
    derivs.registerField(*mInternalDvDx[i]);
    derivs.registerField(*mPairAccelerations[i]);
  }
}

//------------------------------------------------------------------------------
// Initialize the hydro.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {

  // Initialize the grad h corrrections if needed.
  const TableKernel<Dimension>& W = this->kernel();
  const TableKernel<Dimension>& WPi = this->PiKernel();
  if (mGradhCorrection) {
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
    const FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
    FieldList<Dimension, Scalar> omega = state.scalarFields(HydroFieldNames::omegaGradh);
    computeSPHOmegaGradhCorrection(connectivityMap, this->kernel(), position, H, omega);
    for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr) (*boundItr)->applyFieldListGhostBoundary(omega);
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
               WPi);

  // We depend on the caller knowing to finalize the ghost boundaries!
}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
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

  // A few useful constants we'll use in the following loop.
  typedef typename Timing::Time Time;
  const Scalar W0 = W(0.0, 1.0);

  // The connectivity.
  const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  const size_t numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  const FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
  const FieldList<Dimension, Vector> velocity = state.vectorFields(HydroFieldNames::velocity);
  const FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
  const FieldList<Dimension, Scalar> specificThermalEnergy = state.scalarFields(HydroFieldNames::specificThermalEnergy);
  const FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
  const FieldList<Dimension, Scalar> pressure = state.scalarFields(HydroFieldNames::pressure);
  const FieldList<Dimension, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);
  const FieldList<Dimension, Scalar> positionWeight = state.scalarFields(HydroFieldNames::positionWeight);
  const FieldList<Dimension, Scalar> omega = state.scalarFields(HydroFieldNames::omegaGradh);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(positionWeight.size() == numNodeLists);
  CHECK(omega.size() == numNodeLists);

  // Derivative FieldLists.
  FieldList<Dimension, Scalar> rhoSum = derivatives.scalarFields(ReplaceState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  FieldList<Dimension, Vector> DxDt = derivatives.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position);
  FieldList<Dimension, Scalar> DrhoDt = derivatives.scalarFields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::massDensity);
  FieldList<Dimension, Vector> DvDt = derivatives.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);
  FieldList<Dimension, Scalar> DepsDt = derivatives.scalarFields(IncrementState<Dimension, Field<Dimension, Scalar> >::prefix() + HydroFieldNames::specificThermalEnergy);
  FieldList<Dimension, Tensor> DvDx = derivatives.tensorFields(HydroFieldNames::velocityGradient);
  FieldList<Dimension, Tensor> localDvDx = derivatives.tensorFields(HydroFieldNames::internalVelocityGradient);
  FieldList<Dimension, SymTensor> DHDt = derivatives.symTensorFields(IncrementState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  FieldList<Dimension, SymTensor> Hideal = derivatives.symTensorFields(ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  FieldList<Dimension, Scalar> maxViscousPressure = derivatives.scalarFields(HydroFieldNames::maxViscousPressure);
  FieldList<Dimension, vector<Vector> > pairAccelerations = derivatives.vectorVectorFields(HydroFieldNames::pairAccelerations);
  FieldList<Dimension, Scalar> XSPHWeightSum = derivatives.scalarFields(HydroFieldNames::XSPHWeightSum);
  FieldList<Dimension, Vector> XSPHDeltaV = derivatives.vectorFields(HydroFieldNames::XSPHDeltaV);
  FieldList<Dimension, Scalar> weightedNeighborSum = derivatives.scalarFields(HydroFieldNames::weightedNeighborSum);
  FieldList<Dimension, SymTensor> massSecondMoment = derivatives.symTensorFields(HydroFieldNames::massSecondMoment);
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
  CHECK(XSPHWeightSum.size() == numNodeLists);
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
      const Scalar& mi = mass(nodeListi, i);
      const Vector& vi = velocity(nodeListi, i);
      const Scalar& rhoi = massDensity(nodeListi, i);
      const Scalar& epsi = specificThermalEnergy(nodeListi, i);
      const Scalar& Pi = pressure(nodeListi, i);
      const SymTensor& Hi = H(nodeListi, i);
      const Scalar& ci = soundSpeed(nodeListi, i);
      const Scalar& pwi = positionWeight(nodeListi, i);
      const Scalar& omegai = omega(nodeListi, i);
      const Scalar Hdeti = Hi.Determinant();
      const Scalar safeOmegai = omegai/(omegai*omegai + 1.0e-4);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(omegai > 0.0);
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
      Scalar& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      Vector& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      Scalar& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      SymTensor& massSecondMomenti = massSecondMoment(nodeListi, i);
      Scalar& worki = workFieldi(i);

      // Get the connectivity info for this node.
      const vector< vector<int> >& fullConnectivity = connectivityMap.connectivityForNode(&nodeList, i);

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
              const Scalar& mj = mass(nodeListj, j);
              const Vector& vj = velocity(nodeListj, j);
              const Scalar& rhoj = massDensity(nodeListj, j);
              const Scalar& epsj = specificThermalEnergy(nodeListj, j);
              const Scalar& Pj = pressure(nodeListj, j);
              const SymTensor& Hj = H(nodeListj, j);
              const Scalar& cj = soundSpeed(nodeListj, j);
              const Scalar& pwj = positionWeight(nodeListj, j);
              const Scalar& omegaj = omega(nodeListj, j);
              const Scalar Hdetj = Hj.Determinant();
              const Scalar safeOmegaj = omegaj/(omegaj*omegaj + 1.0e-4);
              CHECK(mj > 0.0);
              CHECK(rhoj > 0.0);
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
              Scalar& XSPHWeightSumj = XSPHWeightSum(nodeListj, j);
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
              const Vector Hetai = Hi*etai.unitVector();
              const std::pair<double, double> WWi = W.kernelAndGradValue(etaMagi, Hdeti);
              const Scalar Wi = WWi.first;
              const Scalar gWi = WWi.second;
              const Vector gradWi = gWi*Hetai;

              const Vector Hetaj = Hj*etaj.unitVector();
              const std::pair<double, double> WWj = W.kernelAndGradValue(etaMagj, Hdetj);
              const Scalar Wj = WWj.first;
              const Scalar gWj = WWj.second;
              const Vector gradWj = gWj*Hetaj;

              // Zero'th and second moment of the node distribution -- used for the
              // ideal H calculation.
              const double rij2 = rij.magnitude2();
              const SymTensor thpt = rij.selfdyad()/(rij2 + 1.0e-10) / FastMath::square(Dimension::pownu12(rij2 + 1.0e-10));
              weightedNeighborSumi += fweightij*pwj*std::abs(gWi);
              weightedNeighborSumj += fweightij*pwi*std::abs(gWj);
              massSecondMomenti += fweightij*pwj*gradWi.magnitude2()*thpt;
              massSecondMomentj += fweightij*pwi*gradWj.magnitude2()*thpt;

              // Contribution to the sum density (only if the same material).
              if (nodeListi == nodeListj) {
                rhoSumi += mj*Wi;
                rhoSumj += mi*Wj;
              }

              // Mass density evolution.
              const Vector vij = vi - vj;
              const double deltaDrhoDti = vij.dot(gradWi);
              const double deltaDrhoDtj = vij.dot(gradWj);
              DrhoDti += deltaDrhoDti;
              DrhoDtj += deltaDrhoDtj;

              // Compute the pair-wise artificial viscosity.
              const pair<Tensor, Tensor> QPiij = Q.Piij(nodeListi, i, nodeListj, j,
                                                        ri, etai, vi, rhoi, ci, Hi,
                                                        rj, etaj, vj, rhoj, cj, Hj);
              const Vector Qacci = 0.5*(QPiij.first *gradWi);
              const Vector Qaccj = 0.5*(QPiij.second*gradWj);
//               const Scalar workQi = 0.5*(QPiij.first *vij).dot(gradWi);
//               const Scalar workQj = 0.5*(QPiij.second*vij).dot(gradWj);
              const Scalar workQi = vij.dot(Qacci);
              const Scalar workQj = vij.dot(Qaccj);
              const Scalar Qi = rhoi*rhoi*(QPiij.first. diagonalElements().maxAbsElement());
              const Scalar Qj = rhoj*rhoj*(QPiij.second.diagonalElements().maxAbsElement());
              maxViscousPressurei = max(maxViscousPressurei, Qi);
              maxViscousPressurej = max(maxViscousPressurej, Qj);

              // Determine an effective pressure including a term to fight the tensile instability.
//             const Scalar fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
              const Scalar fij = mEpsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
              const Scalar Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
              const Scalar Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
              const Scalar Peffi = Pi + Ri;
              const Scalar Peffj = Pj + Rj;

              // Acceleration.
              CHECK(rhoi > 0.0);
              CHECK(rhoj > 0.0);
              const double Prhoi = Peffi/(rhoi*rhoi);
              const double Prhoj = Peffj/(rhoj*rhoj);
              const Vector deltaDvDt = Prhoi*safeOmegai*gradWi + Prhoj*safeOmegaj*gradWj + Qacci + Qaccj;
              DvDti -= mj*deltaDvDt;
              DvDtj += mi*deltaDvDt;

              // Specific thermal energy evolution.
              DepsDti += mj*(Prhoi*deltaDrhoDti + workQi);
              DepsDtj += mi*(Prhoj*deltaDrhoDtj + workQj);
              if (mCompatibleEnergyEvolution) {
                if (i < firstGhostNodei) pairAccelerationsi.push_back(-mj*deltaDvDt);
                if (j < firstGhostNodej) pairAccelerationsj.push_back(mi*deltaDvDt);
              }

              // Velocity gradient.
              const Tensor deltaDvDxi = vij.dyad(gradWi);
              const Tensor deltaDvDxj = vij.dyad(gradWj);
              DvDxi -= mj*deltaDvDxi;
              DvDxj -= mi*deltaDvDxj;
              if (nodeListi == nodeListj) {
                localDvDxi -= mj*deltaDvDxi;
                localDvDxj -= mi*deltaDvDxj;
              }

              // Estimate of delta v (for XSPH).
              if (mXSPH and (nodeListi == nodeListj)) {
                const double fXSPH = max(0.0, min(1.0, abs(vij.dot(rij)*safeInv(vij.magnitude()*rij.magnitude()))));
                CHECK(fXSPH >= 0.0 and fXSPH <= 1.0);
                XSPHWeightSumi += fXSPH*mj/rhoj*Wi;
                XSPHWeightSumj += fXSPH*mi/rhoi*Wj;
                XSPHDeltaVi -= fXSPH*mj/rhoj*Wi*vij;
                XSPHDeltaVj += fXSPH*mi/rhoi*Wj*vij;
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

      // Add the self-contribution to density sum.
      rhoSumi += mi*W0*Hdeti;

      // Finish the continuity equation.
      DrhoDti *= mi*safeOmegai;

      // Finish the thermal energy derivative.
      DepsDti *= safeOmegai;

      // Finish the gradient of the velocity.
      CHECK(rhoi > 0.0);
      DvDxi *= safeOmegai/rhoi;
      localDvDxi *= safeOmegai/rhoi;

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (mXSPH) {
        CHECK(XSPHWeightSumi >= 0.0);
        XSPHWeightSumi += Hdeti*mi/rhoi*W0 + 1.0e-30;
        DxDti = vi + XSPHDeltaVi/XSPHWeightSumi;
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
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    FieldList<Dimension, Vector> accelerations = derivs.vectorFields(IncrementState<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::velocity);
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
FVPMTotalHydroBase<Dimension>::
finalize(const typename Dimension::Scalar time,
         const typename Dimension::Scalar dt,
         DataBase<Dimension>& dataBase,
         State<Dimension>& state,
         StateDerivatives<Dimension>& derivs) {

  // Base class finalization.
  GenericHydro<Dimension>::finalize(time, dt, dataBase, state, derivs);

  // Depending on the mass density advancement selected, we may want to replace the 
  // mass density.
  if (densityUpdate() == RigorousSumDensity) {
    const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
    const FieldList<Dimension, Vector> position = state.vectorFields(HydroFieldNames::position);
    const FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
    const FieldList<Dimension, SymTensor> H = state.symTensorFields(HydroFieldNames::H);
    FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
    computeSPHSumMassDensity(connectivityMap, this->kernel(), position, mass, H, massDensity);
  } else if (densityUpdate() == SumDensity) {
    FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
    FieldList<Dimension, Scalar> massDensitySum = derivs.scalarFields(ReplaceState<Dimension, Field<Dimension, Field<Dimension, Scalar> > >::prefix() + 
                                                                      HydroFieldNames::massDensity);
    massDensity.assignFields(massDensitySum);
  }
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& derivs) {

  // Apply boundary conditions to the basic fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.scalarFields(HydroFieldNames::specificThermalEnergy);
  FieldList<Dimension, Vector> velocity = state.vectorFields(HydroFieldNames::velocity);
  FieldList<Dimension, Scalar> pressure = state.scalarFields(HydroFieldNames::pressure);
  FieldList<Dimension, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);
  FieldList<Dimension, Scalar> positionWeight = state.scalarFields(HydroFieldNames::positionWeight);
  FieldList<Dimension, Scalar> omega = state.scalarFields(HydroFieldNames::omegaGradh);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.scalarFields(HydroFieldNames::specificThermalEnergy + "0");

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(positionWeight);
    (*boundaryItr)->applyFieldListGhostBoundary(omega);
    if (compatibleEnergyEvolution()) (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy0);
  }
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.scalarFields(HydroFieldNames::mass);
  FieldList<Dimension, Scalar> massDensity = state.scalarFields(HydroFieldNames::massDensity);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.scalarFields(HydroFieldNames::specificThermalEnergy);
  FieldList<Dimension, Vector> velocity = state.vectorFields(HydroFieldNames::velocity);
  FieldList<Dimension, Scalar> pressure = state.scalarFields(HydroFieldNames::pressure);
  FieldList<Dimension, Scalar> soundSpeed = state.scalarFields(HydroFieldNames::soundSpeed);
  FieldList<Dimension, Scalar> positionWeight = state.scalarFields(HydroFieldNames::positionWeight);
  FieldList<Dimension, Scalar> omega = state.scalarFields(HydroFieldNames::omegaGradh);

  FieldList<Dimension, Scalar> specificThermalEnergy0;
  if (compatibleEnergyEvolution()) specificThermalEnergy0 = state.scalarFields(HydroFieldNames::specificThermalEnergy + "0");

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->enforceFieldListBoundary(positionWeight);
    (*boundaryItr)->enforceFieldListBoundary(omega);
    if (compatibleEnergyEvolution()) (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy0);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  file.write(mPositionWeight, pathName + "/positionWeight");
  file.write(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.write(mHideal, pathName + "/Hideal");
  file.write(mMassDensitySum, pathName + "/massDensitySum");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.write(mOmegaGradh, pathName + "/omegaGradh");
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");

//   this->artificialViscosity().dumpState(file, pathName + "/Q");

}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
FVPMTotalHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  file.read(mPositionWeight, pathName + "/positionWeight");
  file.read(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  file.read(mHideal, pathName + "/Hideal");
  file.read(mMassDensitySum, pathName + "/massDensitySum");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  file.read(mOmegaGradh, pathName + "/omegaGradh");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");

//   this->artificialViscosity().restoreState(file, pathName + "/Q");
}

}

//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class FVPMTotalHydroBase< Dim<1> >;
}
