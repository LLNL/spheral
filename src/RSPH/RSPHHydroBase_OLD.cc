//---------------------------------Spheral++----------------------------------//
// SPHHydroBase -- modified SPHHydro for large density discontinuities
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
//#include "computeSPHSumMassDensity.hh"
//#include "correctSPHSumMassDensity.hh"
//#include "computeSumVoronoiCellMassDensity.hh"
//#include "computeSPHOmegaGradhCorrection.hh"
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
//#include "Hydro/VolumePolicy.hh"
//#include "Hydro/VoronoiMassDensityPolicy.hh"
//#include "Hydro/SumVoronoiMassDensityPolicy.hh"
//#include "Hydro/SpecificThermalEnergyPolicy.hh"
//#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
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
#include "Utilities/Timer.hh"
#include "Mesh/Mesh.hh"
//#include "CRKSPH/volumeSpacing.hh"

#include "RSPH/RSPHHydroBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

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

// Declare timers
extern Timer TIME_RSPH;
extern Timer TIME_RSPHinitializeStartup;
extern Timer TIME_RSPHregister;
extern Timer TIME_RSPHregisterDerivs;
extern Timer TIME_RSPHpreStepInitialize;
extern Timer TIME_RSPHinitialize;
extern Timer TIME_RSPHfinalizeDerivs;
extern Timer TIME_RSPHghostBounds;
extern Timer TIME_RSPHupdateVol;
extern Timer TIME_RSPHenforceBounds;
extern Timer TIME_RSPHevalDerivs;
extern Timer TIME_RSPHevalDerivs_initial;
extern Timer TIME_RSPHevalDerivs_pairs;
extern Timer TIME_RSPHevalDerivs_final;

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
RSPHHydroBase<Dimension>::
RSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
             DataBase<Dimension>& dataBase,
             ArtificialViscosity<Dimension>& Q,
             const TableKernel<Dimension>& W,
             const double cfl,
             const bool useVelocityMagnitudeForDt,
             const bool XSPH,
             const bool correctVelocityGradient,
             const HEvolutionType HUpdate,
             const double epsTensile,
             const double nTensile,
             const Vector& xmin,
             const Vector& xmax):
  //mAlpha(alpha),
  //mDiffusionCoefficient(diffusionCoefficient),
  //mSumDensityNodeLists(sumDensityNodeLists),
  GenericHydro<Dimension>(Q, cfl, useVelocityMagnitudeForDt),
  mKernel(W),
  //mPiKernel(WPi),
  mSmoothingScaleMethod(smoothingScaleMethod),
  //mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  //mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  //mEvolveTotalEnergy(evolveTotalEnergy),
  //mGradhCorrection(gradhCorrection),
  mXSPH(XSPH),
  mCorrectVelocityGradient(correctVelocityGradient),
  //mSumMassDensityOverAllNodeLists(sumMassDensityOverAllNodeLists),
  //mfilter(filter),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mxmin(xmin),
  mxmax(xmax),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  //mOmegaGradh(FieldStorageType::CopyFields),
  //mSpecificThermalEnergy0(FieldStorageType::CopyFields),
  mEntropy(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
  //mMaxViscousPressure(FieldStorageType::CopyFields),
  //mEffViscousPressure(FieldStorageType::CopyFields),
  //mMassDensityCorrection(FieldStorageType::CopyFields),
  //mViscousWork(FieldStorageType::CopyFields),
  //mMassDensitySum(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldStorageType::CopyFields),
  mMassSecondMoment(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  mLocalM(FieldStorageType::CopyFields),
  //mVolume(FieldStorageType::CopyFields),
  mPairAccelerations(),
  mLastDrhoDx(FieldStorageType::CopyFields),
  mLastDcDx(FieldStorageType::CopyFields),
  mLastDpDx(FieldStorageType::CopyFields),
  mLastDvDx(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  //mOmegaGradh = dataBase.newFluidFieldList(1.0, HydroFieldNames::omegaGradh);
  //mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  mEntropy = dataBase.newFluidFieldList(0.0, HydroFieldNames::entropy);
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  //mMaxViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::maxViscousPressure);
  //mEffViscousPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::effectiveViscousPressure);
  //mMassDensityCorrection = dataBase.newFluidFieldList(0.0, HydroFieldNames::massDensityCorrection);
  //mViscousWork = dataBase.newFluidFieldList(0.0, HydroFieldNames::viscousWork);
  //mMassDensitySum = dataBase.newFluidFieldList(0.0, ReplaceFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity);
  mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
  mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
  mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mPairAccelerations.clear();
  mM = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::M_SPHCorrection);
  mLocalM = dataBase.newFluidFieldList(Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection);
  mLastDrhoDx = dataBase.newFluidFieldList(Vector::zero, "densityGradientLastTimeStep");
  mLastDcDx = dataBase.newFluidFieldList(Vector::zero, "soundSpeedGradientLastTimeStep");
  mLastDpDx = dataBase.newFluidFieldList(Vector::zero, "pressureGradientLastTimeStep");
  mLastDvDx = dataBase.newFluidFieldList(Tensor::zero, "velocityGradientLastTimeStep");
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
RSPHHydroBase<Dimension>::
~RSPHHydroBase() {
}


//------------------------------------------------------------------------------
// register our states
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_RSPHregister.start();
  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);

  // FieldLists
  FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();

  // Policies
  PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
  PolicyPointer thermalEnergyPolicy(new IncrementFieldList<Dimension, Scalar>());
  PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position, true));
  PolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
  PolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());

  std::shared_ptr<CompositeFieldListPolicy<Dimension, Scalar> > rhoPolicy(new CompositeFieldListPolicy<Dimension, Scalar>());
  std::shared_ptr<CompositeFieldListPolicy<Dimension, SymTensor> > Hpolicy(new CompositeFieldListPolicy<Dimension, SymTensor>());
  for (typename DataBase<Dimension>::FluidNodeListIterator itr = dataBase.fluidNodeListBegin();
       itr != dataBase.fluidNodeListEnd();
       ++itr) {
    rhoPolicy->push_back(new IncrementBoundedState<Dimension, Scalar>((*itr)->rhoMin(),
                                                                      (*itr)->rhoMax()));
    const Scalar hmaxInv = 1.0/(*itr)->hmax();
    const Scalar hminInv = 1.0/(*itr)->hmin();
    if (HEvolution() == HEvolutionType::IntegrateH) {
      Hpolicy->push_back(new IncrementBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    } else {
      CHECK(HEvolution() == HEvolutionType::IdealH);
      Hpolicy->push_back(new ReplaceBoundedState<Dimension, SymTensor, Scalar>(hmaxInv, hminInv));
    }
  }

  // enroll
  state.enroll(mass);
  state.enroll(massDensity, rhoPolicy);
  state.enroll(Hfield, Hpolicy);
  state.enroll(position, positionPolicy);
  state.enroll(specificThermalEnergy, thermalEnergyPolicy);
  state.enroll(velocity, velocityPolicy);
  state.enroll(mTimeStepMask);
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mSoundSpeed, csPolicy);

  TIME_RSPHregister.stop();
}


//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_RSPHregisterDerivs.start();

  // Create the scratch fields.
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
  //dataBase.resizeFluidFieldList(mMaxViscousPressure, 0.0, HydroFieldNames::maxViscousPressure, false);
  //dataBase.resizeFluidFieldList(mEffViscousPressure, 0.0, HydroFieldNames::effectiveViscousPressure, false);
  //dataBase.resizeFluidFieldList(mMassDensityCorrection, 0.0, HydroFieldNames::massDensityCorrection, false);
  //dataBase.resizeFluidFieldList(mViscousWork, 0.0, HydroFieldNames::viscousWork, false);
  //dataBase.resizeFluidFieldList(mMassDensitySum, 0.0, ReplaceFieldList<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mNormalization, 0.0, HydroFieldNames::normalization, false);
  dataBase.resizeFluidFieldList(mWeightedNeighborSum, 0.0, HydroFieldNames::weightedNeighborSum, false);
  dataBase.resizeFluidFieldList(mMassSecondMoment, SymTensor::zero, HydroFieldNames::massSecondMoment, false);
  dataBase.resizeFluidFieldList(mXSPHWeightSum, 0.0, HydroFieldNames::XSPHWeightSum, false);
  dataBase.resizeFluidFieldList(mXSPHDeltaV, Vector::zero, HydroFieldNames::XSPHDeltaV, false);
  dataBase.resizeFluidFieldList(mDvDt, Vector::zero, HydroFieldNames::hydroAcceleration, false);
  dataBase.resizeFluidFieldList(mDmassDensityDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, false);
  dataBase.resizeFluidFieldList(mDspecificThermalEnergyDt, 0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, false);
  dataBase.resizeFluidFieldList(mDHDt, SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H, false);
  dataBase.resizeFluidFieldList(mDvDx, Tensor::zero, HydroFieldNames::velocityGradient, false);
  dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mM, Tensor::zero, HydroFieldNames::M_SPHCorrection, false);
  dataBase.resizeFluidFieldList(mLocalM, Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection, false);
  derivs.enroll(mHideal);
  //derivs.enroll(mMaxViscousPressure);
  //derivs.enroll(mEffViscousPressure);
  //derivs.enroll(mMassDensityCorrection);
  //derivs.enroll(mViscousWork);
  //derivs.enroll(mMassDensitySum);
  derivs.enroll(mNormalization);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mXSPHWeightSum);
  derivs.enroll(mXSPHDeltaV);

  // Check if someone already registered DxDt.
  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }

  // Check that no-one else is trying to control the hydro vote for DvDt.
  CHECK(not derivs.registered(mDvDt));
  derivs.enroll(mDvDt);

  derivs.enroll(mDmassDensityDt);
  derivs.enroll(mDspecificThermalEnergyDt);
  derivs.enroll(mDHDt);
  derivs.enroll(mDvDx);
  derivs.enroll(mInternalDvDx);
  derivs.enroll(mM);
  derivs.enroll(mLocalM);
  derivs.enrollAny(HydroFieldNames::pairAccelerations, mPairAccelerations);
  TIME_RSPHregisterDerivs.stop();
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {
  TIME_RSPHinitializeStartup.start();
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
  TIME_RSPHinitializeStartup.stop();
}

//------------------------------------------------------------------------------
//  initialize override
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  TIME_RSPHinitialize.start();

  // storage for HLLC reconstruction
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  mLastDvDx = DvDx;
  mLastDvDx.copyFields();

  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr){
            (*boundItr)->applyFieldListGhostBoundary(mLastDvDx);
            //(*boundItr)->applyFieldListGhostBoundary(mLastDpDx);
            //(*boundItr)->applyFieldListGhostBoundary(mLastDrhoDx);
            //(*boundItr)->applyFieldListGhostBoundary(mLastDcDx);
         }
  
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  // We depend on the caller knowing to finalize the ghost boundaries!
  TIME_RSPHinitialize.stop();
}

//------------------------------------------------------------------------------
// R specialized density summmation
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs) {

}

//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
  TIME_RSPHevalDerivs.start();
  TIME_RSPHevalDerivs_initial.start();

  // zero out our HLLC derivs

  //static double totalLoopTime = 0.0;
  // Get the ArtificialViscosity.
  //auto& Q = this->artificialViscosity();

  // The kernels and such.
  const auto& W = this->kernel();
  //const auto& WQ = this->PiKernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);
  //const auto alpha = this->alpha();
  //const auto diffCoeff = this->diffusionCoefficient();
  //const auto evolveTotalEnergy = this->evolveTotalEnergy();
  const auto epsTensile = this->epsilonTensile();
  //const auto compatibleEnergyEvolution = this->compatibleEnergyEvolution();
  const auto XSPH = this->XSPH();
  const auto& smoothingScaleMethod = this->smoothingScaleMethod();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  //const auto omega = state.fields(HydroFieldNames::omegaGradh, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  //CHECK(omega.size() == numNodeLists);

  // Derivative FieldLists.
  //auto  rhoSum = derivatives.fields(ReplaceFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  //auto  normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  //auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivatives.fields("local " + HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  //auto  maxViscousPressure = derivatives.fields(HydroFieldNames::maxViscousPressure, 0.0);
  ///auto  effViscousPressure = derivatives.fields(HydroFieldNames::effectiveViscousPressure, 0.0);
  //auto  viscousWork = derivatives.fields(HydroFieldNames::viscousWork, 0.0);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto  XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  //CHECK(rhoSum.size() == numNodeLists);
  //CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  //CHECK(localDvDx.size() == numNodeLists);
  CHECK(M.size() == numNodeLists);
  //CHECK(localM.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  //CHECK(maxViscousPressure.size() == numNodeLists);
  //CHECK(effViscousPressure.size() == numNodeLists);
  //CHECK(viscousWork.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  pairAccelerations.resize(npairs);

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  TIME_RSPHevalDerivs_initial.stop();

  // Walk all the interacting pairs.
  TIME_RSPHevalDerivs_pairs.start();
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;
    Scalar Wi, gWi, WQi, gWQi, Wj, gWj, WQj, gWQj, Pstar;
    Vector vstar;
    Tensor QPiij, QPiji;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;

    //auto rhoSum_thread = rhoSum.threadCopy(threadStack);
    //auto normalization_thread = normalization.threadCopy(threadStack);
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DrhoDt_thread = DrhoDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    //auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);
    //auto maxViscousPressure_thread = maxViscousPressure.threadCopy(threadStack, ThreadReduction::MAX);
    //auto effViscousPressure_thread = effViscousPressure.threadCopy(threadStack);
    //auto viscousWork_thread = viscousWork.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread = XSPHDeltaV.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      //Vector& DrhoDxi = mLastDrhoDx(nodeListi,i);
      //Vector& DrhoDxj = mLastDrhoDx(nodeListj,j);
      //Vector& DpDxi = mLastDpDx(nodeListi,i);
      //Vector& DpDxj = mLastDpDx(nodeListj,j);
      //Vector& DcDxi = mLastDcDx(nodeListi,i);
      //Vector& DcDxj = mLastDcDx(nodeListj,j);

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi,i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      //const auto& omegai = omega(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      //const auto  safeOmegai = safeInv(omegai, tiny);
      //const auto Ki = max(tiny,rhoi*ci*ci);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      //auto& rhoSumi = rhoSum_thread(nodeListi, i);
      //auto& normi = normalization_thread(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DrhoDti = DrhoDt_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      //auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& Mi = M_thread(nodeListi, i);
      auto& localMi = localM_thread(nodeListi, i);
      //auto& maxViscousPressurei = maxViscousPressure_thread(nodeListi, i);
      //auto& effViscousPressurei = effViscousPressure_thread(nodeListi, i);
      //auto& viscousWorki = viscousWork_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj,j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      //const auto& omegaj = omega(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      //const auto  safeOmegaj = safeInv(omegaj, tiny); 
      //const auto Kj = max(tiny,rhoj*cj*cj);
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      //auto& rhoSumj = rhoSum_thread(nodeListj, j);
      //auto& normj = normalization_thread(nodeListj, j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DrhoDtj = DrhoDt_thread(nodeListj, j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      //auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& Mj = M_thread(nodeListj, j);
      auto& localMj = localM_thread(nodeListj, j);
      //auto& maxViscousPressurej = maxViscousPressure_thread(nodeListj, j);
      //auto& effViscousPressurej = effViscousPressure_thread(nodeListj, j);
      //auto& viscousWorkj = viscousWork_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj, j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij =  (nodeListi == nodeListj);

      // Node displacement.
      const auto rij = ri - rj;
      const auto Hij = 0.5*(Hi+Hj);
      const auto etaij = Hij*rij;
      const auto Hdetij = Hij.Determinant();
      const auto etaMagij = etaij.magnitude();
      CHECK(etaMagij >= 0.0);

      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      std::tie(Wi, gWi) = W.kernelAndGradValue(etaMagi, Hdeti);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      std::tie(Wj, gWj) = W.kernelAndGradValue(etaMagj, Hdetj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      const auto fweightij = sameMatij ? 1.0 : mj*rhoi/(mi*rhoj);
      const auto rij2 = rij.magnitude2();
      const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      weightedNeighborSumi +=     fweightij*std::abs(gWi);
      weightedNeighborSumj += 1.0/fweightij*std::abs(gWj);
      massSecondMomenti +=     fweightij*gradWi.magnitude2()*thpt;
      massSecondMomentj += 1.0/fweightij*gradWj.magnitude2()*thpt;

      // averaged things.
      const auto vij = vi - vj;
      const auto rhoij = 0.5*(rhoi+rhoj); 
      const auto cij = 0.5*(ci+cj);  
      const auto Wij = 0.5*(Wi+Wj); 
      const auto gWij = 0.5*(gWi+gWj);
      const auto gradWij = 0.5*(gradWi+gradWj);
      
      // volumes
      const auto voli = mi/rhoi;
      const auto volj = mj/rhoj;


      // our gradients for HLLC
      //DpDxi += volj*Pi*gradWi;
      //DpDxj += voli*Pj*gradWj;
      //if (sameMatij) {
      //  DrhoDxi += mj*gradWi;
      //  DrhoDxj += mi*gradWj;
      //  DcDxi += volj*cj*gradWi;
      //  DcDxj += voli*ci*gradWj;
      //  localMi -= volj*rij.dyad(gradWi);
      //  localMj -= voli*rij.dyad(gradWj);
      //}

      //std::tie(QPiij, QPiji) = Q.Piij(nodeListi, i, nodeListj, j,
      //                                ri, etaij, vi, rhoij, cij, Hij,  
      //                                rj, etaij, vj, rhoij, cij, Hij); 

      // Determine an effective pressure including a term to fight the tensile instability.
//             const auto fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
      const auto fij = epsTensile*FastMath::pow4(Wij/(Hdeti*WnPerh));
      const auto Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
      const auto Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
      const auto Peffi = Pi + Ri;
      const auto Peffj = Pj + Rj;

      // Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(ci > 0.0);
      CHECK(cj > 0.0);

      const auto isExpanding = vij.dot(rij) > 0.0;

      if (isExpanding){
        vstar = 0.5*(vi+vj);
        Pstar = 0.5*(Pi+Pj);
      }else{
        computeHLLCstate( rij, 
                        nodeListi, nodeListj, i, j,
                        Peffi, Peffj, rhoi, rhoj, vi, vj, ci, cj,
                        vstar, Pstar);
      }
      
      const auto Prho = Pstar/(rhoi*rhoj);
      const auto deltaDvDt = Prho*(gradWij + gradWij);
      DvDti -= mj*deltaDvDt;
      DvDtj += mi*deltaDvDt;

      pairAccelerations[kk] = deltaDvDt;  // Acceleration for i (j anti-symmetric)

      auto deltaDvDxi = (vi-vstar).dyad(gradWij);
      auto deltaDvDxj = (vstar-vj).dyad(gradWij);

      DvDxi -= volj*deltaDvDxi; 
      DvDxj -= voli*deltaDvDxj;

    // Specific thermal energy evolution.
      DepsDti += mj*(Prho*deltaDvDxi.Trace());
      DepsDtj += mi*(Prho*deltaDvDxj.Trace());

      // Linear gradient correction term.
      if(this->mCorrectVelocityGradient){
        Mi -= volj*rij.dyad(gradWij);
        Mj -= voli*rij.dyad(gradWij);
      }
      // Estimate of delta v (for XSPH).
      if (XSPH and (sameMatij)) {
        //const auto wXSPHij = 0.5*(voli*Wij + volj*Wij);
        XSPHWeightSumi += volj*Wij;//wXSPHij;
        XSPHWeightSumj += voli*Wij;//wXSPHij;
        XSPHDeltaVi -= volj*Wij*vij;
        XSPHDeltaVj += voli*Wij*vij;
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region
  TIME_RSPHevalDerivs_pairs.stop();

  // Finish up the derivatives for each point.
  TIME_RSPHevalDerivs_final.start();
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto  hmin = nodeList.hmin();
    const auto  hmax = nodeList.hmax();
    const auto  hminratio = nodeList.hminratio();
    const auto  nPerh = nodeList.nodesPerSmoothingScale();

    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto DpDxi = mLastDpDx(nodeListi,i);
      auto DcDxi = mLastDpDx(nodeListi,i);
      auto DrhoDxi = mLastDpDx(nodeListi,i);

      //auto& rhoSumi = rhoSum(nodeListi, i);
      //auto& normi = normalization(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      //auto& localDvDxi = localDvDx(nodeListi, i);
      //auto& localMi = localM(nodeListi, i);
      auto& Mi = M(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);


      // Finish the gradient of the velocity.
      if (this->mCorrectVelocityGradient and
          std::abs(Mi.Determinant()) > 1.0e-10 and
          numNeighborsi > Dimension::pownu(2)) {
        Mi = Mi.Inverse();
        DvDxi = DvDxi*Mi;
        //DpDxi = Mi*DpDxi;
      } 

      // Finish the gradient of the velocity.
      //if (std::abs(localMi.Determinant()) > 1.0e-10 and
      //    numNeighborsi > Dimension::pownu(2)) {
      //  localMi = localMi.Inverse();
      //  DrhoDxi = localMi*DrhoDxi;
      //  DcDxi = localMi*DcDxi;
      //} 

      // Evaluate the continuity equation.
      DrhoDti -=  rhoi*DvDxi.Trace();

      // If needed finish the total energy derivative.
      //if (evolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      if (XSPH) {
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti = vi + XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      } else {
        DxDti = vi;
      }

      // The H tensor evolution.
      DHDti = smoothingScaleMethod.smoothingScaleDerivative(Hi,
                                                             ri,
                                                             DvDxi,
                                                             hmin,
                                                             hmax,
                                                             hminratio,
                                                             nPerh);
      Hideali = smoothingScaleMethod.newSmoothingScale(Hi,
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
    }
  }
  TIME_RSPHevalDerivs_final.stop();
  TIME_RSPHevalDerivs.stop();
}



template<typename Dimension>
void
RSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {


  // The kernels and such.
  const auto& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists.
  // State FieldLists.
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);

  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);

  // Derivative FieldLists.
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);

  CHECK(M.size() == numNodeLists);


  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  pairAccelerations.resize(npairs);

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  TIME_RSPHevalDerivs_initial.stop();

  // Walk all the interacting pairs.
  TIME_RSPHevalDerivs_pairs.start();
#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;

    auto M_thread = M.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      //Vector& DrhoDxi = mLastDrhoDx(nodeListi,i);
      //Vector& DrhoDxj = mLastDrhoDx(nodeListj,j);
      //Vector& DpDxi = mLastDpDx(nodeListi,i);
      //Vector& DpDxj = mLastDpDx(nodeListj,j);
      //Vector& DcDxi = mLastDcDx(nodeListi,i);
      //Vector& DcDxj = mLastDcDx(nodeListj,j);

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& Mi = M_thread(nodeListi, i);
  
      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);


      auto& Mj = M_thread(nodeListj, j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij =  (nodeListi == nodeListj);

      // Node displacement.
      const auto rij = ri - rj;

      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      const auto gWi = W.gradValue(etaMagi, Hdeti);
      const auto Hetai = Hi*etai.unitVector();
      const auto gradWi = gWi*Hetai;

      const auto gWj = W.gradValue(etaMagj, Hdetj);
      const auto Hetaj = Hj*etaj.unitVector();
      const auto gradWj = gWj*Hetaj;

      // averaged things.
      const auto gradWij = 0.5*(gradWi+gradWj);
      
      // volumes
      const auto voli = mi/rhoi;
      const auto volj = mj/rhoj;

      // Linear gradient correction term.
      if(this->mCorrectVelocityGradient){
        Mij = rij.dyad(gradWij);
        Mi -= volj*Mij;
        Mj -= voli*Mij;
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region

}


//------------------------------------------------------------------------------------
// finish up the derivatives
//------------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_RSPHfinalizeDerivs.start();

  TIME_RSPHfinalizeDerivs.stop();
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {
  TIME_RSPHghostBounds.start();

  // Apply boundary conditions to the basic fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
  }
  TIME_RSPHghostBounds.stop();
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_RSPHenforceBounds.start();

  // Enforce boundary conditions on the fluid state Fields.
  FieldList<Dimension, Scalar> mass = state.fields(HydroFieldNames::mass, 0.0);
  FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  FieldList<Dimension, Scalar> specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  FieldList<Dimension, Vector> velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Scalar> pressure = state.fields(HydroFieldNames::pressure, 0.0);
  FieldList<Dimension, Scalar> soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
  }
  TIME_RSPHenforceBounds.stop();
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  file.write(mTimeStepMask, pathName + "/timeStepMask");
  file.write(mPressure, pathName + "/pressure");
  file.write(mSoundSpeed, pathName + "/soundSpeed");
  //file.write(mVolume, pathName + "/volume");
  //file.write(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  // file.write(mEntropy, pathName + "/entropy");
  file.write(mHideal, pathName + "/Hideal");
  //file.write(mMassDensitySum, pathName + "/massDensitySum");
  //file.write(mNormalization, pathName + "/normalization");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");

  //file.write(mOmegaGradh, pathName + "/omegaGradh");
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
  //file.write(mMaxViscousPressure, pathName + "/maxViscousPressure");
  //file.write(mEffViscousPressure, pathName + "/effectiveViscousPressure");
  //file.write(mMassDensityCorrection, pathName + "/massDensityCorrection");
  //file.write(mViscousWork, pathName + "/viscousWork");
  file.write(mM, pathName + "/M");
  file.write(mLocalM, pathName + "/localM");

//   this->artificialViscosity().dumpState(file, pathName + "/Q");

}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
 
  file.read(mTimeStepMask, pathName + "/timeStepMask");
  file.read(mPressure, pathName + "/pressure");
  file.read(mSoundSpeed, pathName + "/soundSpeed");
  //file.read(mVolume, pathName + "/volume");
  //file.read(mSpecificThermalEnergy0, pathName + "/specificThermalEnergy0");
  // file.read(mEntropy, pathName + "/entropy");
  file.read(mHideal, pathName + "/Hideal");
  //file.read(mMassDensitySum, pathName + "/massDensitySum");
  //file.read(mNormalization, pathName + "/normalization");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  //file.read(mOmegaGradh, pathName + "/omegaGradh");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  //file.read(mMaxViscousPressure, pathName + "/maxViscousPressure");
  file.read(mM, pathName + "/M");
  file.read(mLocalM, pathName + "/localM");

  // For backwards compatibility on change 3597 -- drop in the near future.
  for (auto DvDtPtr: mDvDt) DvDtPtr->name(HydroFieldNames::hydroAcceleration);

//   this->artificialViscosity().restoreState(file, pathName + "/Q");
}

//====================================================================================
//    Pj,rhoj,cj,vj       |         rij          (R)
//    o------------------------------------->o
//     (L)                |->ustar          Pi,rhoi,ci,vi
//====================================================================================
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
computeHLLCstate( const typename Dimension::Vector& rij,
                  int nodeListi,
                  int nodeListj,
                  int i,
                  int j,
                  const typename Dimension::Scalar& Pi,
                  const typename Dimension::Scalar& Pj,
                  const typename Dimension::Scalar& rhoi, 
                  const typename Dimension::Scalar& rhoj,
                  const typename Dimension::Vector& vi,   
                  const typename Dimension::Vector& vj,
                  const typename Dimension::Scalar& ci,   
                  const typename Dimension::Scalar& cj,
                  typename Dimension::Vector& vstar,     
                  typename Dimension::Scalar& Pstar) const {

  const auto DvDxi = mLastDvDx(nodeListi,i);
  const auto DvDxj = mLastDvDx(nodeListj,j);
  const auto DrhoDxi = mLastDrhoDx(nodeListi,i);
  const auto DrhoDxj = mLastDrhoDx(nodeListj,j);
  const auto DpDxi = mLastDpDx(nodeListi,i);
  const auto DpDxj = mLastDpDx(nodeListj,j);
  const auto DcDxi = mLastDcDx(nodeListi,i);
  const auto DcDxj = mLastDcDx(nodeListj,j);

  // normal component
  const auto rhatij = rij.unitVector();
  //const auto isExpanding = (vi-vj).dot(rij) > 0.0;

  // perpendicular component 
  //const auto wij = vi+vj-(ui+uj)*rhatij;                      
  //const auto wi = vi-ui*rhatij;
  //const auto wj = vj-ui*rhatij;
  
  // dumb 1D wave speed 
  //const auto Si = ui + ci;
  //const auto Sj = uj - cj;

  // accoustic impedance
  //const auto Ci = rhoi*(Si-ui);
  //const auto Cj = rhoj*(Sj-uj);

  // get our vanLeer limiter 
  const auto xij = rij/2.0;

  //auto gradi = DrhoDxi.dot(xij);
  //auto gradj = DrhoDxj.dot(xij);
  //auto ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  //auto rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  //auto x = min(ri,rj);

  //auto phi = ( x>0.0 ? 
  //             2.0/(1.0 + x)*2.0*x/(1.0 + x) : 
  //             0.0
  //           );
  
  // project linear soln
  //const auto rhoR = rhoi - phi*DrhoDxi.dot(xij);
  //const auto rhoL = rhoj + phi*DrhoDxj.dot(xij);

  //gradi = DcDxi.dot(xij);
  //gradj = DcDxj.dot(xij);
  //ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  //rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  //x = min(ri,rj);

  //phi = ( x>0.0 ? 
  //        2.0/(1.0 + x)*2.0*x/(1.0 + x) : 
  //        0.0
  //      );
  
  // project linear soln
  //const auto cR = ci - phi*DcDxi.dot(xij);
  //const auto cL = cj + phi*DcDxj.dot(xij);

  //gradi = DpDxi.dot(xij);
  //gradj = DpDxj.dot(xij);
  //ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  //rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  //x = min(ri,rj);

  //phi = ( x>0.0 ? 
  //        2.0/(1.0 + x)*2.0*x/(1.0 + x) : 
  //        0.0
  //      );
  
  // project linear soln
  //const auto pR = Pi - phi*DpDxi.dot(xij);
  //const auto pL = Pj + phi*DpDxj.dot(xij);


  const auto gradi = (DvDxi.dot(xij)).dot(xij);
  const auto gradj = (DvDxj.dot(xij)).dot(xij);
  const auto ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  const auto rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  const auto x = min(ri,rj);

  const auto phi = ( x>0.0 ? 
                     2.0/(1.0 + x)*2.0*x/(1.0 + x) : 
                     0.0
                   );
  
  // project linear soln
  const auto vi1 = vi - phi*DvDxi*xij;
  const auto vj1 = vj + phi*DvDxj*xij;
  const auto ui = vi1.dot(rhatij);
  const auto uj = vj1.dot(rhatij);

    // dont go out of bounds
  //const auto ui = vi.dot(rhatij);
  //const auto uj = vj.dot(rhatij);
  const auto umin = min(ui,uj);
  const auto umax = max(ui,uj);

  // 1D interface velocity 
  const auto ustar = ((ci > 0.0 and cj>0.0) ?
                      min(max((rhoi*ci*ui + rhoj*cj*uj - Pi + Pj )/max((rhoi*ci + rhoj*cj),1e-30),umin),umax) :
                      (ui+uj)/2.0);

  vstar = ustar*rhatij + ((vi+vj) - (ui+uj)*rhatij)/2.0;

  Pstar = rhoj*cj * (uj - ustar) + Pj;

  // that'll do it hopefully my dumb wave speed isn't too dumb
  }

}