//---------------------------------Spheral++----------------------------------//
// RSPHHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/Physics.hh"
#include "SPH/computeSPHSumMassDensity.hh"
//#include "Physics/GenericHydro.hh"
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
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"
//#include "Hydro/EntropyPolicy.hh"
//#include "Mesh/MeshPolicy.hh"
//#include "Mesh/generateMesh.hh"
#include "DataBase/DataBase.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/Timer.hh"
//#include "Mesh/Mesh.hh"
//#include "CRKSPH/volumeSpacing.hh"

#include "FSISpecificThermalEnergyPolicy.hh"
#include "RSPHHydroBase.hh"

#ifdef _OPENMP
#include "omp.h"
#endif

#include <limits.h>
#include <float.h>
#include <algorithm>
#include <fstream>
#include <map>
#include <vector>
#include <sstream>

using std::vector;
using std::string;
using std::pair;
using std::to_string;
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


namespace {
// Provide to_string for Spheral Vector
template<typename Vector>
std::string
vec_to_string(const Vector& vec) {
  std::ostringstream oss;
  oss << vec << std::endl;
  return oss.str();
}
}

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given artificial viscosity and kernels.
//------------------------------------------------------------------------------
template<typename Dimension>
RSPHHydroBase<Dimension>::
RSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
             DataBase<Dimension>& dataBase,
             const TableKernel<Dimension>& W,
             const double cfl,
             const bool useVelocityMagnitudeForDt,
             const bool compatibleEnergyEvolution,
             const bool evolveTotalEnergy,
             const bool XSPH,
             const bool correctVelocityGradient,
             const HEvolutionType HUpdate,
             const double epsTensile,
             const double nTensile,
             const Vector& xmin,
             const Vector& xmax):
  Physics<Dimension>(),
  mCfl(cfl),
  mUseVelocityMagnitudeForDt(useVelocityMagnitudeForDt),
  mMinMasterNeighbor(INT_MAX),
  mMaxMasterNeighbor(0),
  mSumMasterNeighbor(0),
  mMinCoarseNeighbor(INT_MAX),
  mMaxCoarseNeighbor(0),
  mSumCoarseNeighbor(0),
  mMinRefineNeighbor(INT_MAX),
  mMaxRefineNeighbor(0),
  mSumRefineNeighbor(0),
  mMinActualNeighbor(INT_MAX),
  mMaxActualNeighbor(0),
  mSumActualNeighbor(0),
  mNormMasterNeighbor(0),
  mNormCoarseNeighbor(0),
  mNormRefineNeighbor(0),
  mNormActualNeighbor(0),
  mKernel(W),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mCorrectVelocityGradient(correctVelocityGradient),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mxmin(xmin),
  mxmax(xmax),
  mTimeStepMask(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mSpecificThermalEnergy0(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
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
  mPairAccelerations(),
  mDpDx(FieldStorageType::CopyFields),
  mDrhoDx(FieldStorageType::CopyFields),
  mLastDrhoDx(FieldStorageType::CopyFields),
  mLastDpDx(FieldStorageType::CopyFields),
  mLastDvDx(FieldStorageType::CopyFields),
  mLastInternalDvDx(FieldStorageType::CopyFields),
  mRestart(registerWithRestart(*this)) {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mSpecificThermalEnergy0 = dataBase.newFluidFieldList(0.0, HydroFieldNames::specificThermalEnergy + "0");
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
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
  mDrhoDx = dataBase.newFluidFieldList(Vector::zero, "densityGradient");
  mLastDrhoDx = dataBase.newFluidFieldList(Vector::zero, "densityGradientLastTimeStep");
  mDpDx = dataBase.newFluidFieldList(Vector::zero, "pressureGradient");
  mLastDpDx = dataBase.newFluidFieldList(Vector::zero, "pressureGradientLastTimeStep");
  mLastDvDx = dataBase.newFluidFieldList(Tensor::zero, "velocityGradientLastTimeStep");
  mLastInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, "localVelocityGradientLastTimeStep");
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
RSPHHydroBase<Dimension>::
~RSPHHydroBase() {
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
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_RSPHregister.start();

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  // Create the local storage for time step mask, pressure, sound speed, and position weight.
  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);
  dataBase.resizeFluidFieldList(mLastDvDx, Tensor::zero, "velocityGradientLastTimeStep",false);
  dataBase.resizeFluidFieldList(mLastInternalDvDx, Tensor::zero, "localVelocityGradientLastTimeStep",false);
  dataBase.resizeFluidFieldList(mLastDpDx, Vector::zero, "pressureGradientLastTimeStep",false);
  dataBase.resizeFluidFieldList(mLastDrhoDx, Vector::zero, "densityGradientLastTimeStep",false);
  
  //dataBase.resizeFluidFieldList(mSpecificThermalEnergy0, 0.0);

  FieldList<Dimension, Scalar> mass = dataBase.fluidMass();
  FieldList<Dimension, Scalar> massDensity = dataBase.fluidMassDensity();
  FieldList<Dimension, SymTensor> Hfield = dataBase.fluidHfield();
  FieldList<Dimension, Vector> position = dataBase.fluidPosition();
  FieldList<Dimension, Scalar> specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  FieldList<Dimension, Vector> velocity = dataBase.fluidVelocity();
  //FieldList<Dimension, Tensor> DvDx0 = mLastDvDx;
  //FieldList<Dimension, Tensor> DpDx0 = mLastDpDx;

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
  
  PolicyPointer positionPolicy(new IncrementFieldList<Dimension, Vector>());
  PolicyPointer pressurePolicy(new PressurePolicy<Dimension>());
  PolicyPointer csPolicy(new SoundSpeedPolicy<Dimension>());
  
  state.enroll(mTimeStepMask);
  state.enroll(mass);
  //state.enroll(mDpDx);
  state.enroll(massDensity, rhoPolicy);
  state.enroll(Hfield, Hpolicy);
  state.enroll(position, positionPolicy);
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mSoundSpeed, csPolicy);


  // conditional for energy method
  if (mCompatibleEnergyEvolution) {
    
    PolicyPointer thermalEnergyPolicy(new FSISpecificThermalEnergyPolicy<Dimension>(dataBase));
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,
                                                                           HydroFieldNames::specificThermalEnergy,
                                                                           true));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);

  }else if (mEvolveTotalEnergy) {
    PolicyPointer thermalEnergyPolicy(new SpecificFromTotalThermalEnergyPolicy<Dimension>());
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,
                                                                           HydroFieldNames::specificThermalEnergy,
                                                                           true));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);

  } else {
    PolicyPointer thermalEnergyPolicy(new IncrementFieldList<Dimension, Scalar>());
    PolicyPointer velocityPolicy(new IncrementFieldList<Dimension, Vector>(HydroFieldNames::position,
                                                                           true));
    state.enroll(specificThermalEnergy, thermalEnergyPolicy);
    state.enroll(velocity, velocityPolicy);
  }
  
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
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  dataBase.resizeFluidFieldList(mDrhoDx, Vector::zero,"densityGradient", false);
  dataBase.resizeFluidFieldList(mDpDx, Vector::zero,"pressureGradient", false);
  dataBase.resizeFluidFieldList(mHideal, SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H, false);
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

  // Check if someone already registered DxDt.
  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }
  // Check that no-one else is trying to control the hydro vote for DvDt.
  CHECK(not derivs.registered(mDvDt));
  derivs.enroll(mDrhoDx);
  derivs.enroll(mDpDx);
  derivs.enroll(mDvDt);
  derivs.enroll(mHideal);
  derivs.enroll(mNormalization);
  derivs.enroll(mWeightedNeighborSum);
  derivs.enroll(mMassSecondMoment);
  derivs.enroll(mXSPHWeightSum);
  derivs.enroll(mXSPHDeltaV);
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
// RSPH Time step
//------------------------------------------------------------------------------
template<typename Dimension>
typename RSPHHydroBase<Dimension>::TimeStepType
RSPHHydroBase<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const{
     const double tiny = std::numeric_limits<double>::epsilon();

  // Get some useful fluid variables from the DataBase.
  const auto  mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  eps = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  cs = state.fields(HydroFieldNames::soundSpeed, 0.0);
  const auto  DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto& connectivityMap = dataBase.connectivityMap(this->requireGhostConnectivity(),
                                                         this->requireOverlapConnectivity(),
                                                         this->requireIntersectionConnectivity());
  const auto  numNodeLists = connectivityMap.nodeLists().size();

  // Initialize the return value to some impossibly high value.
  auto minDt = make_pair(std::numeric_limits<double>::max(), string());

  // Loop over every fluid node.
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& fluidNodeList = **(dataBase.fluidNodeListBegin() + nodeListi);
    const auto nPerh = fluidNodeList.nodesPerSmoothingScale();
    CHECK(nPerh > 0.0);

    // Walk all the nodes in this FluidNodeList.
    const auto ni = connectivityMap.numNodes(nodeListi);
    // #pragma omp parallel for reduction(MINPAIR:minDt)
#pragma omp parallel
    {
      auto minDt_local = minDt;
#pragma omp for
      for (auto k = 0; k < ni; ++k) {
        const auto i = connectivityMap.ithNode(nodeListi, k);

        // If this node is masked, don't worry about it.
        if (mask(nodeListi, i) == 1) {

          // Get this nodes minimum characteristic smoothing scale.
          CHECK2(H(nodeListi, i).Determinant() >  0.0,
                 "Bad H tensor : " << H(nodeListi, i) << " : " << fluidNodeList.name() << " " << i << " " << fluidNodeList.firstGhostNode());
          const auto& Hi = H(nodeListi, i);
          const Scalar nodeScalei = 1.0/Hi.eigenValues().maxElement()/nPerh;
          // const Scalar nodeScalei = 1.0/Dimension::rootnu(Hi.Determinant());
          //     const Scalar nodeScalei = nodeExtent(nodeListi, i).minElement()/kernelExtent;

          // Sound speed limit.
          const auto csi = cs(nodeListi, i);
          const auto csDt = nodeScalei/(csi + tiny);
          if (csDt < minDt_local.first) {
            minDt_local = make_pair(csDt, ("Sound speed limit: dt = " + to_string(csDt) + "\n" +
                                           "                   cs = " + to_string(cs(nodeListi, i)) + "\n" +
                                           "            nodeScale = " + to_string(nodeScalei) + "\n" +
                                           "             material = " + fluidNodeList.name() + "\n" +
                                           "(nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                           "           @ position = " + vec_to_string(position(nodeListi, i))));
          }

          // Velocity divergence limit.
          const auto divVelocity = DvDx(nodeListi, i).Trace();
          const auto divvDt = 1.0/(std::abs(divVelocity) + tiny);
          if (divvDt < minDt_local.first) {
            minDt_local = make_pair(divvDt, ("Velocity divergence limit: dt = " + to_string(divvDt) + "\n" +
                                             "                 div velocity = " + to_string(divVelocity) + "\n" +
                                             "                     material = " + fluidNodeList.name() + "\n" +
                                             "        (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                             "                   @ position = " + vec_to_string(position(nodeListi, i))));
          }

          // Maximum velocity difference limit.
          //const auto& xi = position(nodeListi, i);
          const auto& vi = velocity(nodeListi, i);
          const auto& fullConnectivity = connectivityMap.connectivityForNode(nodeListi, i);
          for (auto nodeListj = 0u; nodeListj != numNodeLists; ++nodeListj) {
            const auto& connectivity = fullConnectivity[nodeListj];
            for (auto jItr = connectivity.begin();
                 jItr != connectivity.end();
                 ++jItr) {
              const auto  j = *jItr;
              //const auto& xj = position(nodeListj, j);
              const auto& vj = velocity(nodeListj, j);
              const auto& Hj = H(nodeListj, j);
              //const auto  rhoj = rho(nodeListj, j);
              const Scalar nodeScalej = 1.0/Hj.eigenValues().maxElement()/nPerh;
              // const Scalar vij = std::abs((vj - vi).dot((xi - xj).unitVector()));
              const auto  vij = vi - vj;
              const auto  dtVelDiff = std::min(nodeScalei, nodeScalej)*safeInvVar(vij.magnitude(), 1e-30);
              if (dtVelDiff < minDt_local.first) {
                minDt_local = make_pair(dtVelDiff, ("Pairwise velocity difference limit: dt = " + to_string(dtVelDiff) + "\n" + 
                                                    "                              material = " + fluidNodeList.name() + "\n" +
                                                    "                  (nodeListi, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                                    "                  (nodeListj, i, rank) = (" + to_string(nodeListj) + " " + to_string(j) + " " + to_string(Process::getRank()) + ")\n" +
                                                    "                   @ pos(nodeListi, i) = " + vec_to_string(position(nodeListi, i)) + "\n" +
                                                    "                   @ pos(nodeListj, j) = " + vec_to_string(position(nodeListj, j)) + "\n" +
                                                    "                                   vij = " + to_string(vij.magnitude()) + "\n" +
                                                    "                            nodeScalei = " + to_string(nodeScalei) + "\n" +
                                                    "                            nodeScalej = " + to_string(nodeScalej)));
              }

              
            }
          }


          // Total acceleration limit.
          const auto vmagi = vi.magnitude();
          const auto dtAcc = 0.1*std::max(nodeScalei/(vmagi + tiny), vmagi/(DvDt(nodeListi, i).magnitude() + tiny));
          if (dtAcc < minDt_local.first) {
            minDt_local = make_pair(dtAcc, ("Total acceleration limit: dt = " + to_string(dtAcc) + "\n" + 
                                            "              |acceleration| = " + to_string(DvDt(nodeListi, i).magnitude()) + "\n" +
                                            "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                            "                    material = " + fluidNodeList.name() + "\n" +
                                            "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                            "                  @ position = " + vec_to_string(position(nodeListi, i))));
          }

          // If requested, limit against the absolute velocity.
          if (useVelocityMagnitudeForDt()) {
            const auto velDt = nodeScalei/(velocity(nodeListi, i).magnitude() + 1.0e-10);
            if (velDt < minDt_local.first) {
              minDt_local = make_pair(velDt, ("Velocity magnitude limit: dt = " + to_string(velDt) + "\n" +
                                              "                        |vi| = " + to_string(velocity(nodeListi, i).magnitude()) + "\n" +
                                              "                   nodeScale = " + to_string(nodeScalei) + "\n" +
                                              "                    material = " + fluidNodeList.name() + "\n" +
                                              "       (nodeListID, i, rank) = (" + to_string(nodeListi) + " " + to_string(i) + " " + to_string(Process::getRank()) + ")\n" +
                                              "                  @ position = " + vec_to_string(position(nodeListi, i))));
            }
          }
        }
      }

#pragma omp critical
      if (minDt_local.first < minDt.first) minDt = minDt_local;
    }
  }

  // Scale by the cfl safety factor.
  minDt.first *= cfl();
  return minDt;

}
//------------------------------------------------------------------------------
// This method is called once at the beginning of a timestep, after all state registration.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_RSPHpreStepInitialize.start();
  
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  computeSPHSumMassDensity(connectivityMap, this->kernel(), true, position, mass, H, massDensity);
  for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
  for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  
  TIME_RSPHpreStepInitialize.stop();
}

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar /*time*/,
           const typename Dimension::Scalar /*dt*/,
           const DataBase<Dimension>& dataBase,
           State<Dimension>& state,
           StateDerivatives<Dimension>& derivs) {
  TIME_RSPHinitialize.start();

  // store DvDx for use in HLLC reconstruction
  const auto DvDx = derivs.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  const auto localDvDx = derivs.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  const auto DpDx = derivs.fields("pressureGradient", Vector::zero);
  const auto DrhoDx = derivs.fields("densityGradient", Vector::zero);
  const auto  DvDt = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  //auto DvDx0 = state.fields("velocityGradientLastTimeStep",Tensor::zero);
  //DvDx0.copyFields(DvDx);
  //DvDx0 = DvDx;
  //DvDx0.copyFields();
  mLastDrhoDx = DrhoDx;
  mLastDrhoDx.copyFields();
  mLastDvDx = DvDx;
  mLastDvDx.copyFields();
  mLastInternalDvDx = localDvDx;
  mLastInternalDvDx.copyFields();
  mLastDpDx = DvDt;
  mLastDpDx.copyFields();
  mLastDpDx *= rho;
  mLastDpDx *= -1.0;

  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
       boundItr != this->boundaryEnd();
       ++boundItr){
    (*boundItr)->applyFieldListGhostBoundary(mLastDvDx);
    (*boundItr)->applyFieldListGhostBoundary(mLastInternalDvDx);
    (*boundItr)->applyFieldListGhostBoundary(mLastDpDx);
    (*boundItr)->applyFieldListGhostBoundary(mLastDrhoDx);
  }

  
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
 
  TIME_RSPHinitialize.stop();
}


//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
  TIME_RSPHevalDerivs.start();
  TIME_RSPHevalDerivs_initial.start();

  if (this->correctVelocityGradient()){
    this->evaluateSpatialGradients(time,dt,dataBase,state,derivatives);
  }

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
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  const auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);
  //const auto oldDvDx = state.fields("velocityGradientLastTimeStep",Tensor::zero);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(specificThermalEnergy.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(soundSpeed.size() == numNodeLists);
  CHECK(oldDvDx.size() == numNodeLists);

  // Derivative FieldLists.
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto  XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto  DpDx = derivatives.fields("pressureGradient",Vector::zero);
  auto  DrhoDx = derivatives.fields("densityGradient",Vector::zero);
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DpDx.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  //if (mCompatibleEnergyEvolution) pairAccelerations.resize(npairs);
  pairAccelerations.resize(2*npairs);

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
    Scalar Wi, gWi, Wj, gWj, Pstar, rhostari, rhostarj;
    Vector vstar;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DvDt_thread = DvDt.threadCopy(threadStack);
    auto weightedNeighborSum_thread = weightedNeighborSum.threadCopy(threadStack);
    auto massSecondMoment_thread = massSecondMoment.threadCopy(threadStack);
    auto DepsDt_thread = DepsDt.threadCopy(threadStack);
    auto DrhoDt_thread = DrhoDt.threadCopy(threadStack);
    auto DvDx_thread = DvDx.threadCopy(threadStack);
    auto DpDx_thread = DpDx.threadCopy(threadStack);
    auto DrhoDx_thread = DrhoDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& epsi = specificThermalEnergy(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto& ci = soundSpeed(nodeListi, i);
      const auto& oldDvDxi = DvDx(nodeListi,i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& DpDxi = DpDx_thread(nodeListi,i);
      auto& DrhoDxi = DrhoDx_thread(nodeListi,i);
      auto& DrhoDti = DrhoDt_thread(nodeListi, i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      const auto& Mi = M(nodeListi,i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& epsj = specificThermalEnergy(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto& cj = soundSpeed(nodeListj, j);
      const auto& oldDvDxj = DvDx(nodeListj,j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DvDxj = DvDx_thread(nodeListj, j);
      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& DpDxj = DpDx_thread(nodeListj,j);
      auto& DrhoDxj = DrhoDx_thread(nodeListj,j);
      auto& DrhoDtj = DrhoDt_thread(nodeListj, j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);
      const auto& Mj = M(nodeListj,j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij =  (nodeListi == nodeListj);

      // Node displacement.
      const auto rij = ri - rj;
      const auto vij = vi - vj;
      const auto etai = Hi*rij;
      const auto etaj = Hj*rij;
      const auto etaMagi = etai.magnitude();
      const auto etaMagj = etaj.magnitude();
      CHECK(etaMagi >= 0.0);
      CHECK(etaMagj >= 0.0);

      // Symmetrized kernel weight and gradient.
      std::tie(Wi, gWi) = W.kernelAndGradValue(etaMagi, Hdeti);
      const auto Hetai = Hi*etai.unitVector();
      auto gradWi = gWi*Hetai;

      std::tie(Wj, gWj) = W.kernelAndGradValue(etaMagj, Hdetj);
      const auto Hetaj = Hj*etaj.unitVector();
      auto gradWj = gWj*Hetaj;

      const auto Hdetij = (0.5*(Hi+Hj)).Determinant();
      const auto Wij = 0.5*(Wi+Wj);
      const auto gWij = 0.5*(gWi+gWj);
      const auto gradWij = 0.5*(gradWi+gradWj);
      //gradWi = gradWij;
      //gradWj = gradWij;
      //gWi = gWij;
      //gWj = gWij;
      //Wi = Wij;
      //Wj = Wij;
      if (mCorrectVelocityGradient){
        gradWi = Mi.Transpose()*gradWi;
        gradWj = Mj.Transpose()*gradWj;
      }

      // Zero'th and second moment of the node distribution -- used for the
      // ideal H calculation.
      const auto rij2 = rij.magnitude2();
      const auto thpt = rij.selfdyad()*safeInvVar(rij2*rij2*rij2);
      weightedNeighborSumi += std::abs(gWi);
      weightedNeighborSumj += std::abs(gWj);
      massSecondMomenti += gradWi.magnitude2()*thpt;
      massSecondMomentj += gradWj.magnitude2()*thpt;


      // Determine an effective pressure including a term to fight the tensile instability.
      //const auto fij = epsTensile*pow(Wi/(Hdeti*WnPerh), nTensile);
      const auto fij = mEpsTensile*FastMath::pow4(Wi/(Hdeti*WnPerh));
      const auto Ri = fij*(Pi < 0.0 ? -Pi : 0.0);
      const auto Rj = fij*(Pj < 0.0 ? -Pj : 0.0);
      const auto Peffi = Pi + Ri;
      const auto Peffj = Pj + Rj;

      // Acceleration.
      CHECK(rhoi > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(ci > 0.0);
      CHECK(cj > 0.0);
      const auto voli = mi/rhoi;
      const auto volj = mj/rhoj;

      computeHLLCstate( rij, 
                        nodeListi, nodeListj, i, j,
                        Peffi, Peffj, rhoi, rhoj, vi, vj, ci, cj,
                        DpDxi, DpDxj, oldDvDxi, oldDvDxj,
                        vstar, Pstar, rhostari, rhostarj);
      
      // acceleration
      //------------------------------------------------------
      const auto Prhoi = Pstar/(rhostari*rhostarj)*gradWi;
      const auto Prhoj = Pstar/(rhostari*rhostarj)*gradWj;
      const auto deltaDvDt = Prhoi+Prhoj;//Prho*(+gradWij);
      DvDti -= mj*deltaDvDt;
      DvDtj += mi*deltaDvDt;

      pairAccelerations[2*kk] = deltaDvDt; 

      // velocity gradient -> continuity
      //------------------------------------------------------
      const auto deltaDvDxi = (2.0*(vi-vstar).dyad(gradWi));
      const auto deltaDvDxj = (2.0*(vstar-vj).dyad(gradWj));

      DvDxi -= volj*(deltaDvDxi);
      DvDxj -= voli*(deltaDvDxj);
      //if (sameMatij) {
      //  localDvDxi -= volj*deltaDvDxi; 
      //  localDvDxj -= voli*deltaDvDxj;
      //}

      //const auto weighti = deltaDvDxi.Trace();
      //const auto weightj = deltaDvDxj.Trace();
      //const auto deltaV = vij.dot(gradWi+gradWj) - weighti - weightj;
      //const auto wi = abs(weighti)/(abs(weighti)+abs(weightj)+ std::numeric_limits<Scalar>::epsilon());
      
      DrhoDti += rhostari * volj * deltaDvDxi.Trace();
      DrhoDtj += rhostarj * voli * deltaDvDxj.Trace();
      //DrhoDti += rhoi * volj *(weighti+wi*deltaV);
      //DrhoDtj += rhoj * voli *(weightj+(1.0-wi)*deltaV);
      
    // Specific thermal energy evolution.
    //--------------------------------------------------------
      DepsDti += mj*(2.0*Prhoi.dot(vi-vstar));
      DepsDtj += mi*(2.0*Prhoj.dot(vstar-vj));
      pairAccelerations[2*kk+1].x(2.0*Prhoi.dot(vi-vstar)); 
      pairAccelerations[2*kk+1].y(2.0*Prhoj.dot(vstar-vj)); 

      // DrhoDxi -= volj*(rhoi-rhoj)*gradWi;
      // DrhoDxj -= voli*(rhoi-rhoj)*gradWj;
      // DpDxi -= volj*(Pi-Pj)*gradWi;
      // DpDxj -= voli*(Pi-Pj)*gradWj;
      // localDvDxi -= volj*(vi-vj).dyad(gradWi); 
      // localDvDxj -= voli*(vi-vj).dyad(gradWj);


      // diffusions
      //-----------------------------------------------------------
      const auto rhoij = 0.5*(rhoi+rhoj); 
      const auto cij = 0.5*(ci+cj); 
      const auto Hij = 0.5*(Hi+Hj);
      const auto etaij = Hij*rij;
      const auto etaMagij = etaij.magnitude();
      CHECK(etaMagij>0.0);
      if (sameMatij){
        const auto rhoDiffusionCoeff=0.00;
        const auto diffusion =  rhoDiffusionCoeff*(rhoi-rhoj)*cij*etaij.dot(gradWij)/(etaMagij*etaMagij+tiny);
        DrhoDti += volj*diffusion;
        DrhoDtj -= voli*diffusion;
      }

      if (sameMatij){
        const auto epsDiffusionCoeff=0.00;
        const auto diffusion =  epsDiffusionCoeff*(epsi-epsj)*cij*etaij.dot(gradWij)/(rhoij*etaMagij*etaMagij+tiny);
        DepsDti += mj*diffusion;
        DepsDtj -= mi*diffusion;
        pairAccelerations[2*kk+1][0] += mj*diffusion; 
        pairAccelerations[2*kk+1][1] -= mi*diffusion;
      }

      // XSPH
      //-----------------------------------------------------------
      //if (XSPH and sameMatij) {
      //  XSPHWeightSumi += volj*Wi;
      //  XSPHWeightSumj += voli*Wj;
      //  XSPHDeltaVi -= volj*Wi*vij;
      //  XSPHDeltaVj += voli*Wj*vij;
      //}

    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
  } // OpenMP parallel region
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

      auto& normi = normalization(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
     // auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);
      //auto& Mi = M(nodeListi,i);
      // Add the self-contribution to density sum.
      //normi += mi/rhoi*W0*Hdeti;

      // Finish the gradient of the velocity.
      //if (this->mCorrectVelocityGradient and
      //    std::abs(Mi.Determinant()) > 1.0e-10 and
      //    numNeighborsi > Dimension::pownu(2)) {
      //  DvDxi = DvDxi*Mi;
      //} 

      // Evaluate the continuity equation.
      //DrhoDti = - rhoi * DvDxi.Trace();

      

      // If needed finish the total energy derivative.
      if (mEvolveTotalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      DxDti = vi;
      //if (mXSPH) {
      //  XSPHWeightSumi += Hdeti*mi/rhoi*W0;
      //  CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
      //  DxDti += XSPHDeltaVi/max(tiny, XSPHWeightSumi);
      //} 

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
    }
  }
  TIME_RSPHevalDerivs_final.stop();
  TIME_RSPHevalDerivs.stop();
}


//------------------------------------------------------------------------------
// EvalDerivs subroutine for spatial derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
evaluateSpatialGradients(const typename Dimension::Scalar /*time*/,
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
  const auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(velocity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);

  // Derivative FieldLists.
  auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  DpDx = derivatives.fields("pressureGradient",Vector::zero);
  auto  DrhoDx = derivatives.fields("densityGradient",Vector::zero);
  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  localM = derivatives.fields("local"+HydroFieldNames::M_SPHCorrection, Tensor::zero);
  CHECK(M.size() == numNodeLists);
  CHECK(localM.size() == numNodeLists);
  CHECK(DpDx.size() == numNodeLists);
  CHECK(DrhoDx.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

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
    auto DpDx_thread = DpDx.threadCopy(threadStack);
    auto DrhoDx_thread = DrhoDx.threadCopy(threadStack);
    auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto M_thread = M.threadCopy(threadStack);
    auto localM_thread = localM.threadCopy(threadStack);

#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;

      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& vi = velocity(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& DpDxi = DpDx_thread(nodeListi,i);
      auto& DrhoDxi = DrhoDx_thread(nodeListi,i);
      auto& Mi = M_thread(nodeListi, i);
      auto& localMi = localM_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& vj = velocity(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& DpDxj = DpDx_thread(nodeListj,j);
      auto& DrhoDxj = DrhoDx_thread(nodeListj,j);
      auto& Mj = M_thread(nodeListj, j);
      auto& localMj = localM_thread(nodeListj, j);

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
      
      DpDxi -= volj*(Pi-Pj)*gradWi;
      DpDxj -= voli*(Pi-Pj)*gradWj;
      localDvDxi -= volj*(vi-vj).dyad(gradWi); 
      localDvDxj -= voli*(vi-vj).dyad(gradWj);

      // Linear gradient correction term.
      Mi -= volj*rij.dyad(gradWi);
      Mj -= voli*rij.dyad(gradWj);
    
      if (sameMatij){
        DrhoDxi -= volj*(rhoi-rhoj)*gradWi;
        DrhoDxj -= voli*(rhoi-rhoj)*gradWj;
        localMi -= volj*rij.dyad(gradWi);
        localMj -= voli*rij.dyad(gradWj);
      }

    } // loop over pairs

    // Reduce the thread values to the master.
    threadReduceFieldLists<Dimension>(threadStack);

  }   // OpenMP parallel region

  // Finish up the spatial gradient calculation
  for (auto nodeListi = 0u; nodeListi < numNodeLists; ++nodeListi) {
    const auto& nodeList = mass[nodeListi]->nodeList();
    const auto ni = nodeList.numInternalNodes();
#pragma omp parallel for
    for (auto i = 0u; i < ni; ++i) {
        const auto  numNeighborsi = connectivityMap.numNeighborsForNode(nodeListi, i);
        auto& Mi = M(nodeListi, i);
        auto& localMi = localM(nodeListi, i);
        auto& localDvDxi = localDvDx(nodeListi,i);
        auto& DpDxi = DpDx(nodeListi,i);
        auto& DrhoDxi = DrhoDx(nodeListi,i);

        const auto Mdeti = std::abs(Mi.Determinant());
        const auto localMdeti = std::abs(localMi.Determinant());

        const auto enoughNeighbors =  numNeighborsi > Dimension::pownu(2);
        const auto goodM =  (Mdeti > 0.05 and enoughNeighbors);                   
        const auto goodLocalM =  (localMdeti > enoughNeighbors);

        // interp var to blend out M when it gets ill conditioned detM>0.5 away from 1.0
        const auto x = min(1.0, max(0.0, 1.0-2.0*Mdeti));
        const auto localx = min(1.0, max(0.0, 1.0-2.0*localMdeti));

        Mi = (goodM? (1.0-x)*Mi.Inverse() + x*Tensor::one : Tensor::one);
        localMi = (goodLocalM? (1.0-localx)*localMi.Inverse() + localx*Tensor::one : Tensor::one);
        
        if(mCorrectVelocityGradient){
          localDvDxi = localDvDxi*Mi;
          DpDxi = Mi.Transpose()*DpDxi;
          DrhoDxi = localMi.Transpose()*DrhoDxi;
        }
      }
  }

  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr){
           (*boundItr)->applyFieldListGhostBoundary(localM);
           (*boundItr)->applyFieldListGhostBoundary(M);
           (*boundItr)->applyFieldListGhostBoundary(DpDx);
           (*boundItr)->applyFieldListGhostBoundary(DrhoDx);
           (*boundItr)->applyFieldListGhostBoundary(localDvDx);
  }
  
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
 
}
//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
RSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_RSPHfinalizeDerivs.start();
  // If we're using the compatible energy discretization, we need to enforce
  // boundary conditions on the accelerations.
  if (compatibleEnergyEvolution()) {
    auto accelerations = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(accelerations);
    
    for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
         boundaryItr != this->boundaryEnd();
         ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
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
  file.write(mHideal, pathName + "/Hideal");
  file.write(mNormalization, pathName + "/normalization");
  file.write(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.write(mMassSecondMoment, pathName + "/massSecondMoment");
  file.write(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.write(mXSPHDeltaV, pathName + "/XSPHDeltaV");


  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mInternalDvDx, pathName + "/internalDvDx");
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
  file.read(mHideal, pathName + "/Hideal");
  file.read(mNormalization, pathName + "/normalization");
  file.read(mWeightedNeighborSum, pathName + "/weightedNeighborSum");
  file.read(mMassSecondMoment, pathName + "/massSecondMoment");
  file.read(mXSPHWeightSum, pathName + "/XSPHWeightSum");
  file.read(mXSPHDeltaV, pathName + "/XSPHDeltaV");
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mM, pathName + "/M");
  file.read(mLocalM, pathName + "/localM");

  // For backwards compatibility on change 3597 -- drop in the near future.
  for (auto DvDtPtr: mDvDt) DvDtPtr->name(HydroFieldNames::hydroAcceleration);

}


template<typename Dimension>
void
RSPHHydroBase<Dimension>::
pearlLimiter( const typename Dimension::Scalar& Ci,
              const typename Dimension::Scalar& Cj,
              const typename Dimension::Vector& rij,
              const typename Dimension::Vector& vi,   
              const typename Dimension::Vector& vj,
              const typename Dimension::Tensor& DvDxi,
              const typename Dimension::Tensor& DvDxj,
                    typename Dimension::Vector& vstari,   
                    typename Dimension::Vector& vstarj) const{
  const auto xij = rij/2.0;
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
  vstari = vi - phi*DvDxi*xij;
  vstarj = vj + phi*DvDxj*xij;
}

template<typename Dimension>
const typename Dimension::Scalar
RSPHHydroBase<Dimension>::
vanLeerLimiter( const typename Dimension::Vector& xij,
                //const typename Dimension::Vector& vi,   
                //const typename Dimension::Vector& vj,
                const typename Dimension::Tensor& DvDxi,
                const typename Dimension::Tensor& DvDxj) const{
                //      typename Dimension::Vector& vstari,   
                //      typename Dimension::Vector& vstarj) const{
  const auto gradi = (DvDxi.dot(xij)).dot(xij);
  const auto gradj = (DvDxj.dot(xij)).dot(xij);
  const auto ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  const auto rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  const auto x = min(ri,rj);

  const auto phi = ( x>0.0 ? 
                     2.0/(1.0 + x)*2.0*x/(1.0 + x) : 
                     0.0
                   );
  return phi;
  // project linear soln
  //vstari = vi - phi*DvDxi*xij;
  //vstarj = vj + phi*DvDxj*xij;
}

template<typename Dimension>
const typename Dimension::Scalar 
RSPHHydroBase<Dimension>::
vanLeerLimiter( const typename Dimension::Vector& xij,
                //const typename Dimension::Scalar& vi,   
                //const typename Dimension::Scalar& vj,
                const typename Dimension::Vector& DvDxi,
                const typename Dimension::Vector& DvDxj) const{
                //      typename Dimension::Scalar& vstari,   
                //      typename Dimension::Scalar& vstarj) const{
  const auto gradi = DvDxi.dot(xij);
  const auto gradj = DvDxj.dot(xij);
  const auto ri = gradi/(sgn(gradj)*max(1.0e-30, abs(gradj)));
  const auto rj = gradj/(sgn(gradi)*max(1.0e-30, abs(gradi)));
  const auto x = min(ri,rj);

  const auto phi = ( x>0.0 ? 
                     2.0/(1.0 + x)*2.0*x/(1.0 + x) : 
                     0.0
                   );
  return phi;
  // project linear soln
  //vstari = vi - phi*DvDxi.dot(xij);
  //vstarj = vj + phi*DvDxj.dot(xij);
}

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
                  const typename Dimension::Vector& /*DpDxi*/,
                  const typename Dimension::Vector& /*DpDxj*/,
                  const typename Dimension::Tensor& /*DvDxi*/,
                  const typename Dimension::Tensor& /*DvDxj*/,
                  typename Dimension::Vector& vstar,     
                  typename Dimension::Scalar& Pstar,
                  typename Dimension::Scalar& rhostariOut,
                  typename Dimension::Scalar& rhostarjOut) const {

  // default to average values for bad sound speed
  vstar = (vi+vj)/2.0;
  Pstar = (Pi+Pj)/2.0;
  rhostariOut = rhoi;
  rhostarjOut = rhoj;

  const auto goodC = ( ci > 0.0 and cj > 0.0);

  if (goodC){
    // linear interpolate R and L state
    const auto DvDxi = mLastDvDx(nodeListi,i);
    const auto DvDxj = mLastDvDx(nodeListj,j);
    const auto localDvDxi = mLastDvDx(nodeListi,i);
    const auto localDvDxj = mLastDvDx(nodeListj,j);
    const auto DpDxi = mLastDpDx(nodeListi,i);
    const auto DpDxj = mLastDpDx(nodeListj,j);
    const auto DrhoDxi = mLastDrhoDx(nodeListi,i);
    const auto DrhoDxj = mLastDrhoDx(nodeListj,j);

    const auto xij = rij/2.0;
    const auto phirho = this->vanLeerLimiter( xij, DrhoDxi, DrhoDxj);
    const auto phiv   = this->vanLeerLimiter( xij, localDvDxi, localDvDxj);
    const auto phip   = this->vanLeerLimiter( xij, DpDxi, DpDxj);
    const auto phi = max(phiv,phip);
    const auto vstari = vi - phiv*DvDxi.dot(xij);
    const auto vstarj = vj + phiv*DvDxj.dot(xij);
    const auto pstari = Pi - phip*DpDxi.dot(xij);
    const auto pstarj = Pj + phip*DpDxj.dot(xij);
    const auto rhostari = rhoi;// - phirho*DrhoDxi.dot(xij);
    const auto rhostarj = rhoj;// + phirho*DrhoDxj.dot(xij);

    // get our components
    const auto rhatij = rij.unitVector();

    const auto ui = vstari.dot(rhatij);
    const auto uj = vstarj.dot(rhatij);
    const auto wi = vstari - ui*rhatij;
    const auto wj = vstarj - uj*rhatij;

    // get our wave speed (default acoustic)
    auto Si =   rhostari*ci;
    auto Sj = - rhostarj*cj;
    if (true){ // Davis Einfeldt 1988 min/max treatment
      Si = rhostari*(max(uj+cj,ui+ci)-ui);
      Sj = rhostarj*(min(ui-ci,uj-cj)-uj);
    }else if(false){ // Davis Einfeldt 1988 min/max treatment modified for multi mat
      Si = rhostari*(max(uj-ui,0.0)+ci);
      Sj = rhostarj*(min(ui-uj,0.0)-cj);
    }else if(false){ // HLLE Roe avg +- wavespeed estimate
      const auto sRhoi = sqrt(rhostari);
      const auto sRhoj = sqrt(rhostarj);
      const auto denom = safeInv(sRhoi+sRhoj);
      const auto eta = 0.5*sRhoi*sRhoj*denom*denom;
      const auto d2 = (sRhoi * ci*ci + sRhoj * cj*cj)*denom + eta * (ui-uj)*(ui-uj);
      const auto d = sqrt(d2);

      const auto utilde = (sRhoi * ui + sRhoj * uj)*denom;

      Si = rhostari*((utilde + d)-ui);
      Sj = rhostarj*((utilde - d)-uj);
    }

    // construct our ARS
    const auto denom = safeInv(Si - Sj);
    const auto ustar = (Si*ui - Sj*uj - pstari + pstarj )*denom;
    const auto wstar = (Si*wi - Sj*wj)*denom;
    Pstar = Sj * (ustar-uj) + pstarj;
    vstar = ustar*rhatij + wstar;

    rhostariOut = rhoi;//rhostari*(ui-Si)/max(ustar-Si,1e-30);
    rhostarjOut = rhoj;//rhostarj*(uj-Sj)/max(ustar-Sj,1e-30);
  }
}

}


// Trying out RSPH to reduce the number of neighbors required to resolve a sedov problem. 
// gets negative STE when zero pressure. slightly greater than one good to go.
// XSPH might help sedov w/ ARSPH ? 
// different limiters? 
// other wave speed estimates? HLLE min max type deal
