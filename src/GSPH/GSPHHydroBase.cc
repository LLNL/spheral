//---------------------------------Spheral++----------------------------------//
// GSPHHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include "Physics/Physics.hh"
#include "FileIO/FileIO.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "SPH/computeSPHSumMassDensity.hh"

#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceFieldList.hh"
#include "DataBase/IncrementBoundedFieldList.hh"
#include "DataBase/ReplaceBoundedFieldList.hh"
#include "DataBase/IncrementBoundedState.hh"
#include "DataBase/ReplaceBoundedState.hh"
#include "DataBase/CompositeFieldListPolicy.hh"

#include "Hydro/HydroFieldNames.hh"
#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
#include "Hydro/CompatibleDifferenceSpecificThermalEnergyPolicy.hh"
#include "Hydro/PositionPolicy.hh"
#include "Hydro/PressurePolicy.hh"
#include "Hydro/SoundSpeedPolicy.hh"

#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/timingUtilities.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/globalBoundingVolumes.hh"
#include "Utilities/Timer.hh"

#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/GSPHHydroBase.hh"
#include "GSPH/computeSumVolume.hh"
#include "GSPH/computeSPHVolume.hh"
#include "GSPH/computeMFMDensity.hh"
#include "GSPH/RiemannSolvers/RiemannSolverBase.hh"

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
extern Timer TIME_GSPH;
extern Timer TIME_GSPHinitializeStartup;
extern Timer TIME_GSPHregister;
extern Timer TIME_GSPHregisterDerivs;
extern Timer TIME_GSPHpreStepInitialize;
extern Timer TIME_GSPHinitialize;
extern Timer TIME_GSPHfinalizeDerivs;
extern Timer TIME_GSPHghostBounds;
extern Timer TIME_GSPHenforceBounds;
extern Timer TIME_GSPHevalDerivs;


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
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GSPHHydroBase<Dimension>::
GSPHHydroBase(const SmoothingScaleBase<Dimension>& smoothingScaleMethod,
             DataBase<Dimension>& dataBase,
             RiemannSolverBase<Dimension>& riemannSolver,
             const TableKernel<Dimension>& W,
             const Scalar epsDiffusionCoeff,
             const double cfl,
             const bool useVelocityMagnitudeForDt,
             const bool compatibleEnergyEvolution,
             const bool evolveTotalEnergy,
             const bool XSPH,
             const bool correctVelocityGradient,
             const MassDensityType densityUpdate,
             const HEvolutionType HUpdate,
             const double epsTensile,
             const double nTensile,
             const Vector& xmin,
             const Vector& xmax):
  Physics<Dimension>(),
  mRestart(registerWithRestart(*this)),
  mRiemannSolver(riemannSolver),
  mKernel(W),
  mSmoothingScaleMethod(smoothingScaleMethod),
  mDensityUpdate(densityUpdate),
  mHEvolution(HUpdate),
  mCompatibleEnergyEvolution(compatibleEnergyEvolution),
  mEvolveTotalEnergy(evolveTotalEnergy),
  mXSPH(XSPH),
  mCorrectVelocityGradient(correctVelocityGradient),
  mUseVelocityMagnitudeForDt(useVelocityMagnitudeForDt),
  mEpsTensile(epsTensile),
  mnTensile(nTensile),
  mSpecificThermalEnergyDiffusionCoefficient(epsDiffusionCoeff),
  mCfl(cfl),
  mxmin(xmin),
  mxmax(xmax),
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
  mTimeStepMask(FieldStorageType::CopyFields),
  mVolume(FieldStorageType::CopyFields),
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldStorageType::CopyFields),
  mMassSecondMoment(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mLocalM(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mInternalDvDx(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  mDvDxRaw(FieldStorageType::CopyFields),
  mDpDx(FieldStorageType::CopyFields),
  mDpDxRaw(FieldStorageType::CopyFields),
  mPairAccelerations(),
  mPairDepsDt() {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mVolume = dataBase.newFluidFieldList(0.0, HydroFieldNames::volume);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
  mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
  mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mLocalM = dataBase.newFluidFieldList(Tensor::zero, "local" + HydroFieldNames::M_SPHCorrection);
  mM = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::M_SPHCorrection);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H);
  mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  mDvDxRaw = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient+"RAW");
  mDpDx = dataBase.newFluidFieldList(Vector::zero, GSPHFieldNames::pressureGradient);
  mDpDxRaw = dataBase.newFluidFieldList(Vector::zero, GSPHFieldNames::pressureGradient+"RAW");
  mPairAccelerations.clear();
  mPairDepsDt.clear();

  auto& DpDx = mRiemannSolver.DpDx();
  auto& DvDx = mRiemannSolver.DvDx();
  DpDx = dataBase.newFluidFieldList(Vector::zero,GSPHFieldNames::RiemannPressureGradient);
  DvDx = dataBase.newFluidFieldList(Tensor::zero,GSPHFieldNames::RiemannVelocityGradient);
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
GSPHHydroBase<Dimension>::
~GSPHHydroBase() {
}

//------------------------------------------------------------------------------
// On problem start up, we need to initialize our internal data.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  TIME_GSPHinitializeStartup.start();
  dataBase.fluidPressure(mPressure);
  dataBase.fluidSoundSpeed(mSoundSpeed);
  TIME_GSPHinitializeStartup.stop();
}

//------------------------------------------------------------------------------
// Register the state we need/are going to evolve.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {
  TIME_GSPHregister.start();

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  VERIFY2(not (mCompatibleEnergyEvolution and mEvolveTotalEnergy),
          "SPH error : you cannot simultaneously use both compatibleEnergyEvolution and evolveTotalEnergy");

  dataBase.resizeFluidFieldList(mTimeStepMask, 1, HydroFieldNames::timeStepMask);

  auto mass = dataBase.fluidMass();
  auto massDensity = dataBase.fluidMassDensity();
  auto Hfield = dataBase.fluidHfield();
  auto position = dataBase.fluidPosition();
  auto specificThermalEnergy = dataBase.fluidSpecificThermalEnergy();
  auto velocity = dataBase.fluidVelocity();
  
  auto& DpDx = mRiemannSolver.DpDx();
  auto& DvDx = mRiemannSolver.DvDx();

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

  // normal state variables
  state.enroll(mTimeStepMask);
  state.enroll(mVolume);
  state.enroll(mass);
  state.enroll(massDensity, rhoPolicy);
  state.enroll(Hfield, Hpolicy);
  state.enroll(position, positionPolicy);
  state.enroll(mPressure, pressurePolicy);
  state.enroll(mSoundSpeed, csPolicy);
  state.enroll(DpDx);
  state.enroll(DvDx);

  // conditional for energy method
  if (mCompatibleEnergyEvolution) {
    
    PolicyPointer thermalEnergyPolicy(new CompatibleDifferenceSpecificThermalEnergyPolicy<Dimension>(dataBase));
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
  
  TIME_GSPHregister.stop();
}

//------------------------------------------------------------------------------
// Register the state derivative fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
registerDerivatives(DataBase<Dimension>& dataBase,
                    StateDerivatives<Dimension>& derivs) {
  TIME_GSPHregisterDerivs.start();

  // Create the scratch fields.
  // Note we deliberately do not zero out the derivatives here!  This is because the previous step
  // info here may be used by other algorithms (like the CheapSynchronousRK2 integrator or
  dataBase.resizeFluidFieldList(mDpDx, Vector::zero, GSPHFieldNames::pressureGradient, false);
  dataBase.resizeFluidFieldList(mDpDxRaw, Vector::zero, GSPHFieldNames::pressureGradient+"RAW", false);
  dataBase.resizeFluidFieldList(mDvDxRaw, Tensor::zero, HydroFieldNames::velocityGradient+"RAW", false);
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
  dataBase.resizeFluidFieldList(mLocalM, Tensor::zero, "local" + HydroFieldNames::M_SPHCorrection, false);

  // Check if someone already registered DxDt.
  if (not derivs.registered(mDxDt)) {
    dataBase.resizeFluidFieldList(mDxDt, Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, false);
    derivs.enroll(mDxDt);
  }
  // Check that no-one else is trying to control the hydro vote for DvDt.
  CHECK(not derivs.registered(mDvDt));
  derivs.enroll(mDpDx);
  derivs.enroll(mDpDxRaw);
  derivs.enroll(mDvDxRaw);
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
  derivs.enrollAny(HydroFieldNames::pairWork, mPairDepsDt);
  TIME_GSPHregisterDerivs.stop();
}

//------------------------------------------------------------------------------
// GSPH Time step
//------------------------------------------------------------------------------
template<typename Dimension>
typename GSPHHydroBase<Dimension>::TimeStepType
GSPHHydroBase<Dimension>::
dt(const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar /*currentTime*/) const{
  
  const auto tiny = std::numeric_limits<Scalar>::epsilon();

  // Get some useful fluid variables from the DataBase.
  const auto  mask = state.fields(HydroFieldNames::timeStepMask, 1);
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  const auto  rho = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto  eps = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
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
GSPHHydroBase<Dimension>::
preStepInitialize(const DataBase<Dimension>& dataBase, 
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_GSPHpreStepInitialize.start();
  if(mDensityUpdate == MassDensityType::RigorousSumDensity){
    const auto& connectivityMap = dataBase.connectivityMap();
    const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
    const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
    const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
    auto        massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
    computeSPHSumMassDensity(connectivityMap, this->kernel(), true, position, mass, H, massDensity);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    for (auto boundaryItr = this->boundaryBegin(); boundaryItr < this->boundaryEnd(); ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();
  }
  TIME_GSPHpreStepInitialize.stop();
}

//------------------------------------------------------------------------------
// Initialize the hydro before calling evaluateDerivatives
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
initialize(const typename Dimension::Scalar time,
           const typename Dimension::Scalar dt,
           const DataBase<Dimension>& dataBase,
                 State<Dimension>& state,
                 StateDerivatives<Dimension>& derivs) {
  TIME_GSPHinitialize.start();

  auto& riemannSolver = this->riemannSolver();
  const auto& W = this->kernel();

  riemannSolver.initialize(dataBase, 
                           state,
                           derivs,
                           this->boundaryBegin(),
                           this->boundaryEnd(),
                           time, 
                           dt,
                           W);

  // plop into an intialize volume function
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto  position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto  H = state.fields(HydroFieldNames::H, SymTensor::zero);
  const auto  mass = state.fields(HydroFieldNames::mass, 0.0);
        auto  massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
        auto  volume = state.fields(HydroFieldNames::volume, 0.0);
  if(true){
    computeSumVolume(connectivityMap,W,position,H,volume);
    computeMFMDensity(mass,volume,massDensity);
  }else{
    computeSPHVolume(mass,massDensity,volume);
  }
  
  for (auto boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr){
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
  }
  for (auto boundaryItr = this->boundaryBegin(); 
       boundaryItr < this->boundaryEnd(); 
       ++boundaryItr) (*boundaryItr)->finalizeGhostBoundary();

  TIME_GSPHinitialize.stop();
  
}

//------------------------------------------------------------------------------
// Finalize the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
finalizeDerivatives(const typename Dimension::Scalar /*time*/,
                    const typename Dimension::Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& /*state*/,
                    StateDerivatives<Dimension>& derivs) const {
  TIME_GSPHfinalizeDerivs.start();
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
  TIME_GSPHfinalizeDerivs.stop();
}

//------------------------------------------------------------------------------
// Apply the ghost boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {
  TIME_GSPHghostBounds.start();

  // Apply boundary conditions to the basic fluid state Fields.
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  // our store vars in the riemann solver
  auto DpDx = state.fields(GSPHFieldNames::RiemannPressureGradient,Vector::zero); 
  auto DvDx = state.fields(GSPHFieldNames::RiemannVelocityGradient,Tensor::zero); 

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyFieldListGhostBoundary(volume);
    (*boundaryItr)->applyFieldListGhostBoundary(mass);
    (*boundaryItr)->applyFieldListGhostBoundary(massDensity);
    (*boundaryItr)->applyFieldListGhostBoundary(specificThermalEnergy);
    (*boundaryItr)->applyFieldListGhostBoundary(velocity);
    (*boundaryItr)->applyFieldListGhostBoundary(pressure);
    (*boundaryItr)->applyFieldListGhostBoundary(soundSpeed);
    (*boundaryItr)->applyFieldListGhostBoundary(DpDx);
    (*boundaryItr)->applyFieldListGhostBoundary(DvDx);

  }
  TIME_GSPHghostBounds.stop();
}

//------------------------------------------------------------------------------
// Enforce the boundary conditions for hydro state fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {
  TIME_GSPHenforceBounds.start();

  // Enforce boundary conditions on the fluid state Fields.
  auto volume = state.fields(HydroFieldNames::volume, 0.0);
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  // our store vars in the riemann solver
  auto DpDx = state.fields(GSPHFieldNames::RiemannPressureGradient,Vector::zero); 
  auto DvDx = state.fields(GSPHFieldNames::RiemannVelocityGradient,Tensor::zero); 


  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceFieldListBoundary(volume);
    (*boundaryItr)->enforceFieldListBoundary(mass);
    (*boundaryItr)->enforceFieldListBoundary(massDensity);
    (*boundaryItr)->enforceFieldListBoundary(specificThermalEnergy);
    (*boundaryItr)->enforceFieldListBoundary(velocity);
    (*boundaryItr)->enforceFieldListBoundary(pressure);
    (*boundaryItr)->enforceFieldListBoundary(soundSpeed);
    (*boundaryItr)->enforceFieldListBoundary(DpDx);
    (*boundaryItr)->enforceFieldListBoundary(DvDx);
  }
  TIME_GSPHenforceBounds.stop();
}


//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
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

  // time derivs
  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");

  // spatial derivs
  file.write(mM, pathName + "/M");
  file.write(mDvDx, pathName + "/DvDx");
  file.write(mDpDx, pathName + "/DpDx");
  file.write(mDvDxRaw, pathName + "/DvDxRaw");
  file.write(mDpDxRaw, pathName + "/DpDxRaw");
  

}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
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

  // time derivs
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");

  // spatial derivs
  file.read(mM, pathName + "/M");
  file.read(mDvDx, pathName + "/DvDx");
  file.read(mDpDx, pathName + "/DpDx");
  file.read(mDvDxRaw, pathName + "/DvDxRaw");
  file.read(mDpDxRaw, pathName + "/DpDxRaw");

}

}
