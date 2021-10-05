//---------------------------------Spheral++----------------------------------//
// GSPHHydroBase -- The SPH/ASPH hydrodynamic package for Spheral++.
//
// Created by JMO, Mon Jul 19 22:11:09 PDT 2010
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "NodeList/SmoothingScaleBase.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Physics/Physics.hh"
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

#include "Hydro/SpecificFromTotalThermalEnergyPolicy.hh"
#include "Hydro/SpecificThermalEnergyPolicy.hh"
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

#include "FSISPH/FSIFieldNames.hh"
#include "GSPH/GSPHFieldNames.hh"
#include "GSPH/GSPHHydroBase.hh"
#include "GSPH/GSPHSpecificThermalEnergyPolicy.hh"
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
extern Timer TIME_GSPHupdateVol;
extern Timer TIME_GSPHenforceBounds;
extern Timer TIME_GSPHevalDerivs;
extern Timer TIME_GSPHevalDerivs_initial;
extern Timer TIME_GSPHevalDerivs_pairs;
extern Timer TIME_GSPHevalDerivs_final;


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
  mPressure(FieldStorageType::CopyFields),
  mSoundSpeed(FieldStorageType::CopyFields),
  mHideal(FieldStorageType::CopyFields),
  mNormalization(FieldStorageType::CopyFields),
  mWeightedNeighborSum(FieldStorageType::CopyFields),
  mMassSecondMoment(FieldStorageType::CopyFields),
  mXSPHWeightSum(FieldStorageType::CopyFields),
  mXSPHDeltaV(FieldStorageType::CopyFields),
  mM(FieldStorageType::CopyFields),
  //mLocalM(FieldStorageType::CopyFields),
  mDxDt(FieldStorageType::CopyFields),
  mDvDt(FieldStorageType::CopyFields),
  mDmassDensityDt(FieldStorageType::CopyFields),
  mDspecificThermalEnergyDt(FieldStorageType::CopyFields),
  mDHDt(FieldStorageType::CopyFields),
  mDvDx(FieldStorageType::CopyFields),
  //mInternalDvDx(FieldStorageType::CopyFields),
  mDvDxRaw(FieldStorageType::CopyFields),
  mDpDx(FieldStorageType::CopyFields),
  mDpDxRaw(FieldStorageType::CopyFields),
  mPairAccelerations(),
  mPairDepsDt() {

  // Create storage for our internal state.
  mTimeStepMask = dataBase.newFluidFieldList(int(0), HydroFieldNames::timeStepMask);
  mPressure = dataBase.newFluidFieldList(0.0, HydroFieldNames::pressure);
  mSoundSpeed = dataBase.newFluidFieldList(0.0, HydroFieldNames::soundSpeed);
  mHideal = dataBase.newFluidFieldList(SymTensor::zero, ReplaceBoundedState<Dimension, Field<Dimension, SymTensor> >::prefix() + HydroFieldNames::H);
  mNormalization = dataBase.newFluidFieldList(0.0, HydroFieldNames::normalization);
  mWeightedNeighborSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::weightedNeighborSum);
  mMassSecondMoment = dataBase.newFluidFieldList(SymTensor::zero, HydroFieldNames::massSecondMoment);
  mXSPHWeightSum = dataBase.newFluidFieldList(0.0, HydroFieldNames::XSPHWeightSum);
  mXSPHDeltaV = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::XSPHDeltaV);
  mM = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::M_SPHCorrection);
  //mLocalM = dataBase.newFluidFieldList(Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection);
  mDxDt = dataBase.newFluidFieldList(Vector::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position);
  mDvDt = dataBase.newFluidFieldList(Vector::zero, HydroFieldNames::hydroAcceleration);
  mDmassDensityDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity);
  mDspecificThermalEnergyDt = dataBase.newFluidFieldList(0.0, IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy);
  mDHDt = dataBase.newFluidFieldList(SymTensor::zero, IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::H);
  mDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient);
  //mInternalDvDx = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::internalVelocityGradient);
  mDvDxRaw = dataBase.newFluidFieldList(Tensor::zero, HydroFieldNames::velocityGradient+"RAW");
  mDpDx = dataBase.newFluidFieldList(Vector::zero, GSPHFieldNames::pressureGradient);
  mDpDxRaw = dataBase.newFluidFieldList(Vector::zero, GSPHFieldNames::pressureGradient+"RAW");
  mPairAccelerations.clear();
  mPairDepsDt.clear();

  auto& DpDx = mRiemannSolver.DpDx();
  auto& DvDx = mRiemannSolver.DvDx();
  DpDx = dataBase.newFluidFieldList(Vector::zero,GSPHFieldNames::pressureGradient+"Riemann");
  DvDx = dataBase.newFluidFieldList(Tensor::zero,HydroFieldNames::velocityGradient+"Riemann");
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
    
    PolicyPointer thermalEnergyPolicy(new GSPHSpecificThermalEnergyPolicy<Dimension>(dataBase));
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
  //dataBase.resizeFluidFieldList(mInternalDvDx, Tensor::zero, HydroFieldNames::internalVelocityGradient, false);
  dataBase.resizeFluidFieldList(mM, Tensor::zero, HydroFieldNames::M_SPHCorrection, false);
  //dataBase.resizeFluidFieldList(mLocalM, Tensor::zero, "local " + HydroFieldNames::M_SPHCorrection, false);

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
  //derivs.enroll(mInternalDvDx);
  derivs.enroll(mM);
  //derivs.enroll(mLocalM);
  derivs.enrollAny(HydroFieldNames::pairAccelerations, mPairAccelerations);
  derivs.enrollAny(FSIFieldNames::pairDepsDt, mPairDepsDt);
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
  const TableKernel<Dimension>& W = this->kernel();

  riemannSolver.initialize(dataBase, 
                           state,
                           derivs,
                           this->boundaryBegin(),
                           this->boundaryEnd(),
                           time, 
                           dt,
                           W);
  
  TIME_GSPHinitialize.stop();
  
}


//------------------------------------------------------------------------------
// Determine the principle derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
evaluateDerivatives(const typename Dimension::Scalar time,
                    const typename Dimension::Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivatives) const {
  TIME_GSPHevalDerivs.start();
  TIME_GSPHevalDerivs_initial.start();

  if (this->correctVelocityGradient()) this->evaluateSpatialGradients(time,dt,dataBase,state,derivatives);

  // The kernels and such.
  const auto& W = this->kernel();
  const auto& riemannSolver = this->riemannSolver();

  // A few useful constants we'll use in the following loop.
  const auto tiny = std::numeric_limits<Scalar>::epsilon();
  const auto W0 = W(0.0, 1.0);
  const auto xsph = this->XSPH();
  const auto epsDiffusionCoeff = this->specificThermalEnergyDiffusionCoefficient();
  const auto compatibleEnergy = this->compatibleEnergyEvolution();
  const auto totalEnergy = this->evolveTotalEnergy();

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
  const auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  auto  normalization = derivatives.fields(HydroFieldNames::normalization, 0.0);
  auto  DxDt = derivatives.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::position, Vector::zero);
  auto  DrhoDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::massDensity, 0.0);
  auto  DvDt = derivatives.fields(HydroFieldNames::hydroAcceleration, Vector::zero);
  auto  DepsDt = derivatives.fields(IncrementFieldList<Dimension, Scalar>::prefix() + HydroFieldNames::specificThermalEnergy, 0.0);
  auto  DvDx = derivatives.fields(HydroFieldNames::velocityGradient, Tensor::zero);
  //auto  localDvDx = derivatives.fields(HydroFieldNames::internalVelocityGradient, Tensor::zero);
  auto  DHDt = derivatives.fields(IncrementFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto  Hideal = derivatives.fields(ReplaceBoundedFieldList<Dimension, SymTensor>::prefix() + HydroFieldNames::H, SymTensor::zero);
  auto& pairAccelerations = derivatives.getAny(HydroFieldNames::pairAccelerations, vector<Vector>());
  auto& pairDepsDt = derivatives.getAny(FSIFieldNames::pairDepsDt, vector<Scalar>());
  auto  XSPHWeightSum = derivatives.fields(HydroFieldNames::XSPHWeightSum, 0.0);
  auto  XSPHDeltaV = derivatives.fields(HydroFieldNames::XSPHDeltaV, Vector::zero);
  auto  weightedNeighborSum = derivatives.fields(HydroFieldNames::weightedNeighborSum, 0.0);
  auto  massSecondMoment = derivatives.fields(HydroFieldNames::massSecondMoment, SymTensor::zero);
  auto  DpDx = derivatives.fields(GSPHFieldNames::pressureGradient,Vector::zero);
  auto  DpDxRaw = derivatives.fields(GSPHFieldNames::pressureGradient+"RAW",Vector::zero);
  auto  DvDxRaw = derivatives.fields(HydroFieldNames::velocityGradient+"RAW",Tensor::zero);
  
  CHECK(normalization.size() == numNodeLists);
  CHECK(DxDt.size() == numNodeLists);
  CHECK(DrhoDt.size() == numNodeLists);
  CHECK(DvDt.size() == numNodeLists);
  CHECK(DepsDt.size() == numNodeLists);
  CHECK(DvDx.size() == numNodeLists);
  //CHECK(localDvDx.size() == numNodeLists);
  CHECK(DHDt.size() == numNodeLists);
  CHECK(Hideal.size() == numNodeLists);
  CHECK(XSPHWeightSum.size() == numNodeLists);
  CHECK(XSPHDeltaV.size() == numNodeLists);
  CHECK(weightedNeighborSum.size() == numNodeLists);
  CHECK(massSecondMoment.size() == numNodeLists);
  CHECK(DpDx.size() == numNodeLists);
  //CHECK(DrhoDx.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // Size up the pair-wise accelerations before we start.
  if (compatibleEnergy){
    pairAccelerations.resize(npairs);
    pairDepsDt.resize(2*npairs);
  }

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);
  TIME_GSPHevalDerivs_initial.stop();

  // Walk all the interacting pairs.
  TIME_GSPHevalDerivs_pairs.start();
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
    auto DpDxRaw_thread = DpDxRaw.threadCopy(threadStack);
    auto DvDxRaw_thread = DvDxRaw.threadCopy(threadStack);
    //auto DrhoDx_thread = DrhoDx.threadCopy(threadStack);
    //auto localDvDx_thread = localDvDx.threadCopy(threadStack);
    auto XSPHWeightSum_thread = XSPHWeightSum.threadCopy(threadStack);
    auto XSPHDeltaV_thread =  XSPHDeltaV.threadCopy(threadStack);
    auto normalization_thread = normalization.threadCopy(threadStack);
    
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
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& normi = normalization_thread(nodeListi,i);
      auto& DepsDti = DepsDt_thread(nodeListi, i);
      auto& DvDti = DvDt_thread(nodeListi, i);
      auto& DpDxi = DpDx_thread(nodeListi,i);
      auto& DpDxRawi = DpDxRaw_thread(nodeListi,i);
      auto& DvDxRawi = DvDxRaw_thread(nodeListi,i);
      auto& DvDxi = DvDx_thread(nodeListi, i);
      //auto& localDvDxi = localDvDx_thread(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum_thread(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment_thread(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum_thread(nodeListi,i);
      auto& XSPHDeltaVi = XSPHDeltaV_thread(nodeListi,i);
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
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& normj = normalization_thread(nodeListj,j);
      auto& DvDtj = DvDt_thread(nodeListj, j);
      auto& DepsDtj = DepsDt_thread(nodeListj, j);
      auto& DpDxj = DpDx_thread(nodeListj,j);
      auto& DpDxRawj = DpDxRaw_thread(nodeListj,j);
      auto& DvDxRawj = DvDxRaw_thread(nodeListj,j);
      auto& DvDxj = DvDx_thread(nodeListj, j);
      //auto& localDvDxj = localDvDx_thread(nodeListj, j);
      auto& weightedNeighborSumj = weightedNeighborSum_thread(nodeListj, j);
      auto& massSecondMomentj = massSecondMoment_thread(nodeListj, j);
      auto& XSPHWeightSumj = XSPHWeightSum_thread(nodeListj,j);
      auto& XSPHDeltaVj = XSPHDeltaV_thread(nodeListj,j);
      const auto& Mj = M(nodeListj,j);

      // Flag if this is a contiguous material pair or not.
      const bool sameMatij =  (nodeListi == nodeListj);

      // Node displacement.
      const auto rij = ri - rj;
      const auto rhatij =rij.unitVector();
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

      // const auto Hij = 0.5*(Hi+Hj);
      // const auto Hdetij = Hij.Determinant();
      // const auto etaij = Hij*rij;
      // const auto etaMagij = etaij.magnitude();

      // std::tie(Wij, gWij) = W.kernelAndGradValue(etaMagij, Hdetij);
      // const auto Hetaij = Hij*etaij.unitVector();
      // const auto gradWij = gWij*Hetaij;

      // const auto gradWij = 0.5 * (gradWi + gradWj);
      // const auto gWij = 0.5*(gWi+gWj);
      // const auto Wij = 0.5*(Wi+Wj);
      // Wi = Wij;
      // Wj = Wij;
      // gWi = gWij;
      // gWj = gWij;
      // gradWi = gradWij;
      // gradWj = gradWij;

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

      riemannSolver.interfaceState(i,         j, 
                                   nodeListi, nodeListj, 
                                   ri,        rj, 
                                   rhoi,      rhoj, 
                                   ci,        cj, 
                                   Peffi,     Peffj, 
                                   vi,        vj,    
                                   Pstar,
                                   vstar,
                                   rhostari,
                                   rhostarj);
      
      normi += volj*Wi;
      normj += voli*Wj;

      // acceleration
      //------------------------------------------------------
      const auto Prhoi = Pstar/(rhoi*rhoj)*gradWi;
      const auto Prhoj = Pstar/(rhoi*rhoj)*gradWj;
      const auto deltaDvDt = Prhoi+Prhoj;
      DvDti -= mj*deltaDvDt;
      DvDtj += mi*deltaDvDt;

      // gradients
      //------------------------------------------------------
      const auto deltaDvDxi = 2.0*(vi-vstar).dyad(gradWi);
      const auto deltaDvDxj = 2.0*(vstar-vj).dyad(gradWj);

      // based on riemann soln
      DvDxi -= volj*(deltaDvDxi);
      DvDxj -= voli*(deltaDvDxj);

      DpDxi -= 2.0*volj*(Peffi-Pstar)*gradWi;
      DpDxj -= 2.0*voli*(Pstar-Peffj)*gradWj;

      // based on nodal values
      DvDxRawi -= volj*(vi-vj).dyad(gradWi);
      DvDxRawj -= voli*(vi-vj).dyad(gradWj);

      DpDxRawi -= volj*(Pi-Pj)*gradWi;
      DpDxRawj -= voli*(Pi-Pj)*gradWj;
      
      // specific thermal energy evolution.
      //--------------------------------------------------------
      const auto deltaDepsDti = Prhoi.dot(vi-vstar);
      const auto deltaDepsDtj = Prhoj.dot(vstar-vj);
      DepsDti += mj*(deltaDepsDti);
      DepsDtj += mi*(deltaDepsDtj);
     
      if(compatibleEnergy){
        pairAccelerations[kk] = deltaDvDt; 
        pairDepsDt[2*kk]   = deltaDepsDti; 
        pairDepsDt[2*kk+1] = deltaDepsDtj; 
      }

      // diffusion
      //-----------------------------------------------------------
      if (sameMatij and epsDiffusionCoeff>tiny){
        const auto rhoij = 0.5*(rhoi+rhoj); 
        const auto cij = 0.5*(ci+cj); 
        const auto gradWij = 0.5*(gradWi+gradWj);
        const auto vMagij = vij.dot(rhatij);
        const auto cijEff = max(min(cij + vMagij, cij),0.0);
        const auto diffusion =  epsDiffusionCoeff*(epsi-epsj)*cijEff*rhatij.dot(gradWij)/(rhoij+tiny);
        pairDepsDt[2*kk]   += diffusion; 
        pairDepsDt[2*kk+1] -= diffusion;
      }

      // XSPH
      //-----------------------------------------------------------
      if (xsph) {
        XSPHWeightSumi += volj*Wi;
        XSPHWeightSumj += voli*Wj;
        XSPHDeltaVi -= volj*Wi*(vi-vstar);
        XSPHDeltaVj -= voli*Wj*(vj-vstar);
      }

    } // loop over pairs
    threadReduceFieldLists<Dimension>(threadStack);
  } // OpenMP parallel region
  TIME_GSPHevalDerivs_pairs.stop();


  // Finish up the derivatives for each point.
  TIME_GSPHevalDerivs_final.start();
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
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& normi = normalization(nodeListi, i);
      auto& DxDti = DxDt(nodeListi, i);
      auto& DrhoDti = DrhoDt(nodeListi, i);
      auto& DvDti = DvDt(nodeListi, i);
      auto& DepsDti = DepsDt(nodeListi, i);
      auto& DvDxi = DvDx(nodeListi, i);
      //auto& localDvDxi = localDvDx(nodeListi, i);
      auto& DHDti = DHDt(nodeListi, i);
      auto& Hideali = Hideal(nodeListi, i);
      auto& XSPHWeightSumi = XSPHWeightSum(nodeListi, i);
      auto& XSPHDeltaVi = XSPHDeltaV(nodeListi, i);
      auto& weightedNeighborSumi = weightedNeighborSum(nodeListi, i);
      auto& massSecondMomenti = massSecondMoment(nodeListi, i);

      normi += mi/rhoi*Hdeti*W0;
      //normi /= Hdeti;
      //normi *= 3.1415;

      DrhoDti = - rhoi * DvDxi.Trace() ;

      // If needed finish the total energy derivative.
      if (totalEnergy) DepsDti = mi*(vi.dot(DvDti) + DepsDti);

      // Complete the moments of the node distribution for use in the ideal H calculation.
      weightedNeighborSumi = Dimension::rootnu(max(0.0, weightedNeighborSumi/Hdeti));
      massSecondMomenti /= Hdeti*Hdeti;

      // Determine the position evolution, based on whether we're doing XSPH or not.
      DxDti = vi;
      if (xsph){
        XSPHWeightSumi += Hdeti*mi/rhoi*W0;
        CHECK2(XSPHWeightSumi != 0.0, i << " " << XSPHWeightSumi);
        DxDti += 0.25*XSPHDeltaVi/max(tiny, XSPHWeightSumi);
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
    }
  }
  TIME_GSPHevalDerivs_final.stop();
  TIME_GSPHevalDerivs.stop();
}


//------------------------------------------------------------------------------
// EvalDerivs subroutine for spatial derivs
//------------------------------------------------------------------------------
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
evaluateSpatialGradients(const typename Dimension::Scalar /*time*/,
                         const typename Dimension::Scalar /*dt*/,
                         const DataBase<Dimension>& dataBase,
                         const State<Dimension>& state,
                               StateDerivatives<Dimension>& derivatives) const {


  // The kernels and such.
  const auto& W = this->kernel();

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  // Get the state and derivative FieldLists. 
  const auto mass = state.fields(HydroFieldNames::mass, 0.0);
  const auto position = state.fields(HydroFieldNames::position, Vector::zero);
  const auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const auto H = state.fields(HydroFieldNames::H, SymTensor::zero);
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);

  auto  M = derivatives.fields(HydroFieldNames::M_SPHCorrection, Tensor::zero);
  CHECK(M.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

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

      // const auto Hij = 0.5*(Hi+Hj);
      // const auto Hdetij = Hij.Determinant();
      // const auto etaij = Hij*rij;
      // const auto etaMagij = etaij.magnitude();
      // const auto gWij = W.gradValue(etaMagij, Hdetij);
      // const auto Hetaij = Hij*etaij.unitVector();
      // const auto gradWij = gWij*Hetaij;
      //const auto gradWij = 0.5*(gradWi+gradWj);

      // Linear gradient correction term.
      Mi -= mj/rhoj*rij.dyad(gradWi);
      Mj -= mi/rhoi*rij.dyad(gradWj);
      
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

      const auto Mdeti = std::abs(Mi.Determinant());

      const auto enoughNeighbors =  numNeighborsi > Dimension::pownu(2);
      const auto goodM =  (Mdeti > 1.0e-10 and enoughNeighbors);                   

      // vanLeer limiter for M correction
      // const auto x = min(1.0, max(0.0, 8.0*Mdeti/((1.0+Mdeti)*(1.0+Mdeti)) ));
      // stretched vanleer in a crop top
      // const auto x = min(1.0, max(0.0, 40.0*Mdeti/max((1.0+Mdeti)*(1.0+Mdeti),1.0e-30) ));
      Mi = ( goodM ? Mi.Inverse() : Tensor::one);//( goodM ? x*Mi.Inverse() + (1.0-x)*Tensor::one  : Tensor::one );
    }
    
  }
  
  for (ConstBoundaryIterator boundItr = this->boundaryBegin();
         boundItr != this->boundaryEnd();
         ++boundItr){
           (*boundItr)->applyFieldListGhostBoundary(M);
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
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  // our store vars in the riemann solver
  auto DpDx = state.fields(GSPHFieldNames::pressureGradient+"Riemann",Vector::zero); 
  auto DvDx = state.fields(HydroFieldNames::velocityGradient+"Riemann",Tensor::zero); 

  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
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
  auto mass = state.fields(HydroFieldNames::mass, 0.0);
  auto massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  auto specificThermalEnergy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  auto velocity = state.fields(HydroFieldNames::velocity, Vector::zero);
  auto pressure = state.fields(HydroFieldNames::pressure, 0.0);
  auto soundSpeed = state.fields(HydroFieldNames::soundSpeed, 0.0);

  // our store vars in the riemann solver
  auto DpDx = state.fields(GSPHFieldNames::pressureGradient+"Riemann",Vector::zero); 
  auto DvDx = state.fields(HydroFieldNames::velocityGradient+"Riemann",Tensor::zero); 


  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
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


  file.write(mDxDt, pathName + "/DxDt");
  file.write(mDvDt, pathName + "/DvDt");
  file.write(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.write(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.write(mDHDt, pathName + "/DHDt");
  file.write(mDvDx, pathName + "/DvDx");
  //file.write(mInternalDvDx, pathName + "/internalDvDx");
  file.write(mM, pathName + "/M");
  //file.write(mLocalM, pathName + "/localM");

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
  file.read(mDxDt, pathName + "/DxDt");
  file.read(mDvDt, pathName + "/DvDt");
  file.read(mDmassDensityDt, pathName + "/DmassDensityDt");
  file.read(mDspecificThermalEnergyDt, pathName + "/DspecificThermalEnergyDt");
  file.read(mDHDt, pathName + "/DHDt");
  file.read(mDvDx, pathName + "/DvDx");
  //file.read(mInternalDvDx, pathName + "/internalDvDx");
  file.read(mM, pathName + "/M");
  //file.read(mLocalM, pathName + "/localM");


}


/*
template<typename Dimension>
void
GSPHHydroBase<Dimension>::
computeInitialGradients(Dimension::State<Dimension> state) const {

  // The kernels and such.
  const auto& W = this->kernel();

  // A few useful constants we'll use in the following loop.
  const double tiny = 1.0e-30;
  const Scalar W0 = W(0.0, 1.0);

  // The connectivity.
  const auto& connectivityMap = dataBase.connectivityMap();
  const auto& nodeLists = connectivityMap.nodeLists();
  const auto numNodeLists = nodeLists.size();

  auto mass = dataBase.fluidMass();
  auto massDensity = dataBase.fluidMassDensity();
  const auto H = dataBase.fluidHfield();
  const auto position = dataBase.fluidPosition();
  auto pressure = mPressure;
  CHECK(mass.size() == numNodeLists);
  CHECK(position.size() == numNodeLists);
  CHECK(massDensity.size() == numNodeLists);
  CHECK(pressure.size() == numNodeLists);
  CHECK(H.size() == numNodeLists);

  // The set of interacting node pairs.
  const auto& pairs = connectivityMap.nodePairList();
  const auto  npairs = pairs.size();

  // The scale for the tensile correction.
  const auto& nodeList = mass[0]->nodeList();
  const auto  nPerh = nodeList.nodesPerSmoothingScale();
  const auto  WnPerh = W(1.0/nPerh, 1.0);

#pragma omp parallel
  {
    // Thread private scratch variables
    int i, j, nodeListi, nodeListj;

    typename SpheralThreads<Dimension>::FieldListStack threadStack;
    auto DpDx_thread = mDpDx.threadCopy(threadStack);
    auto M_thread = mM.threadCopy(threadStack);


#pragma omp for
    for (auto kk = 0u; kk < npairs; ++kk) {
      i = pairs[kk].i_node;
      j = pairs[kk].j_node;
      nodeListi = pairs[kk].i_list;
      nodeListj = pairs[kk].j_list;
      
      // Get the state for node i.
      const auto& ri = position(nodeListi, i);
      const auto& mi = mass(nodeListi, i);
      const auto& Pi = pressure(nodeListi, i);
      const auto& rhoi = massDensity(nodeListi, i);
      const auto& Hi = H(nodeListi, i);
      const auto  Hdeti = Hi.Determinant();
      CHECK(mi > 0.0);
      CHECK(rhoi > 0.0);
      CHECK(Hdeti > 0.0);

      auto& DpDxi = DpDx_thread(nodeListi,i);
      auto& Mi = M_thread(nodeListi, i);

      // Get the state for node j
      const auto& rj = position(nodeListj, j);
      const auto& mj = mass(nodeListj, j);
      const auto& Pj = pressure(nodeListj, j);
      const auto& rhoj = massDensity(nodeListj, j);
      const auto& Hj = H(nodeListj, j);
      const auto  Hdetj = Hj.Determinant();
      CHECK(mj > 0.0);
      CHECK(rhoj > 0.0);
      CHECK(Hdetj > 0.0);

      auto& DpDxj = DpDx_thread(nodeListj,j);
      auto& Mj = M_thread(nodeListj, j);

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

      // volumes
      const auto voli = mi/rhoi;
      const auto volj = mj/rhoj;
      
      DpDxi -= volj*(Pi-Pj)*gradWi;
      DpDxj -= voli*(Pi-Pj)*gradWj;

      // Linear gradient correction term.
      Mi -= volj*rij.dyad(gradWi);
      Mj -= voli*rij.dyad(gradWj);
      
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
      auto& Mi = mM(nodeListi, i);
      auto& DpDxi = mDpDx(nodeListi,i);
      auto& lastDpDxi = mLastDpDx(nodeListi,i);

      const auto Mdeti = std::abs(Mi.Determinant());

      const auto enoughNeighbors =  numNeighborsi > Dimension::pownu(2);
      const auto goodM =  (Mdeti > 0.05 and enoughNeighbors);                   

      // interp var to blend out M when it gets ill conditioned detM>0.5 away from 1.0
      const auto x = min(1.0, max(0.0, 1.0-2.0*Mdeti));
      Mi = (goodM? (1.0-x)*Mi.Inverse() + x*Tensor::one : Tensor::one);

      if(mCorrectVelocityGradient){
        DpDxi = Mi.Transpose()*DpDxi;
      }

      lastDpDxi = DpDxi;
    }
    
  }
*/
}


// Trying out GSPH to reduce the number of neighbors required to resolve a sedov problem. 
// gets negative STE when zero pressure. slightly greater than one good to go.
// XSPH might help sedov w/ AGSPH ? 
// different limiters? 
// other wave speed estimates? HLLE min max type deal
