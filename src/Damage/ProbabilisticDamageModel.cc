//---------------------------------Spheral++----------------------------------//
// ProbabilisticDamageModel
// A damage model based on Weibull statistics that uses volume based
// probabilities per node to decide when damage starts to accrue.  Should
// generate similar results to the classic Benz-Asphaug (Grady-Kipp) model
// without generating explicit flaws.  Also appropriate for use with varying
// resolution materials.
//
// Created by JMO, Tue Apr 13 15:58:08 PDT 2021
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "ProbabilisticDamageModel.hh"
#include "TensorStrainPolicy.hh"
#include "ProbabilisticDamagePolicy.hh"
#include "YoungsModulusPolicy.hh"
#include "LongitudinalSoundSpeedPolicy.hh"
#include "DamageGradientPolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/ReplaceState.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/Neighbor.hh"
#include "Utilities/mortonOrderIndices.hh"
#include "Utilities/allReduce.hh"

#include <boost/functional/hash.hpp>  // hash_combine

#include <string>
#include <vector>
#include <algorithm>
#include <limits>
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
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ProbabilisticDamageModel<Dimension>::
ProbabilisticDamageModel(SolidNodeList<Dimension>& nodeList,
                         const TableKernel<Dimension>& W,
                         const double kWeibull,
                         const double mWeibull,
                         const size_t seed,
                         const size_t minFlawsPerNode,
                         const double crackGrowthMultiplier,
                         const double volumeMultiplier,
                         const DamageCouplingAlgorithm damageCouplingAlgorithm,
                         const TensorStrainAlgorithm strainAlgorithm,
                         const bool damageInCompression,
                         const double criticalDamageThreshold,
                         const Field<Dimension, int>& mask):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, damageCouplingAlgorithm),
  mStrainAlgorithm(strainAlgorithm),
  mDamageInCompression(damageInCompression),
  mkWeibull(kWeibull),
  mmWeibull(mWeibull),
  mVolumeMultiplier(volumeMultiplier),
  mVmin(std::numeric_limits<double>::max()),
  mVmax(std::numeric_limits<double>::min()),
  mCriticalDamageThreshold(criticalDamageThreshold),
  mSeed(seed),
  mMinFlawsPerNode(minFlawsPerNode),
  mNumFlaws(SolidFieldNames::numFlaws, nodeList),
  mMask(mask),
  mMinFlaw(SolidFieldNames::minFlaw, nodeList),
  mMaxFlaw(SolidFieldNames::maxFlaw, nodeList),
  mInitialVolume(SolidFieldNames::initialVolume, nodeList),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mDdamageDt(ProbabilisticDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList),
  mRandomGenerator(SolidFieldNames::randomGenerator, nodeList) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ProbabilisticDamageModel<Dimension>::
~ProbabilisticDamageModel() {
}

//------------------------------------------------------------------------------
// initializeProblemStartup
//
// After all initial state has been initialize (node positions, masses, etc),
// but before we try to run any physics cycles.  This is when we initialize a
// lot of our state for the damage model.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
initializeProblemStartup(DataBase<Dimension>& dataBase) {

  // How many points are actually being damaged?
  const auto nused_local = mMask.sumElements();
  const auto nused_global = allReduce(nused_local, MPI_SUM, Communicator::communicator());

  // Compute the Morton-ordering for hashing with the global seed to seed each
  // point-wise random number generator.
  typedef KeyTraits::Key Key;
  const FieldList<Dimension, Key> keyList = mortonOrderIndices(dataBase);
  const auto& nodes = this->nodeList();
  CHECK(keyList.fieldForNodeList(nodes) < keyList.end());
  const Field<Dimension, Key>& keys = **(keyList.fieldForNodeList(nodes));
  
  // Compute the initial volumes and random seeds for each node.  We hash the
  // morton index of each point with the global seed value to create unique
  // but reproducible seeds for each points random number generator.
  const auto& mass = nodes.mass();
  const auto& rho = nodes.massDensity();
  const auto  nlocal = nodes.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    if (mMask(i) == 1) {
      CHECK(mass(i) > 0.0 and rho(i) > 0.0);
      mInitialVolume(i) = mass(i)/rho(i) * mVolumeMultiplier;
      mVmin = std::min(mVmin, mInitialVolume(i));
      mVmax = std::max(mVmax, mInitialVolume(i));
      Key seedi = mSeed;
      boost::hash_combine(seedi, keys(i));
      mRandomGenerator(i).seed(seedi);      // starting out generating in [0,1)
      mRandomGenerator(i)();                // Recommended to discard first value in sequence
    }
  }
  mVmin = allReduce(mVmin, MPI_MIN, Communicator::communicator());
  mVmax = allReduce(mVmax, MPI_MAX, Communicator::communicator());

  // Compute the maximum strain we expect for the minimum volume.
  const auto epsMax2m = mMinFlawsPerNode/(mkWeibull*mVmin);  // epsmax ** m

  // Based on this compute the maximum number of flaws any node will have.  We'll use this to
  // spin the random number generator without extra communiction.
  const auto maxFlawsPerNode = int(std::max(1, int(mkWeibull*mVmax*epsMax2m + 0.5)));

  // Generate initial realizations of the flaw population for each point, of which we capture the min/max range
  // for each point.
  const auto mInv = 1.0/mmWeibull;
  size_t minNumFlaws = std::numeric_limits<size_t>::max();
  size_t maxNumFlaws = 0u;
  size_t totalNumFlaws = 0u;
  auto epsMin = std::numeric_limits<double>::max();
  auto epsMax = std::numeric_limits<double>::min();
  auto sumFlaws = 0.0;
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    if (mMask(i) == 1) {
      mMinFlaw(i) = std::numeric_limits<double>::max();
      mNumFlaws(i) = std::max(1, std::min(maxFlawsPerNode, int(mkWeibull*mInitialVolume(i)*epsMax2m + 0.5)));
      const auto Ai = mNumFlaws(i)/(mkWeibull*mInitialVolume(i));
      CHECK(Ai > 0.0);
      auto sumFlawsi = 0.0;
      for (auto j = 0; j < mNumFlaws(i); ++j) {
        const auto flaw = pow(Ai * mRandomGenerator(i)(), mInv);
        mMinFlaw(i) = std::min(mMinFlaw(i), flaw);
        mMaxFlaw(i) = std::max(mMinFlaw(i), flaw);
        sumFlawsi += flaw;
      }

      // Now reset the random number generator so we'll only generate flaws in the allowed range
      // for this point.
      mRandomGenerator(i).range(pow(mMinFlaw(i), mmWeibull)/Ai, pow(mMaxFlaw(i), mmWeibull)/Ai);
      CHECK(mRandomGenerator(i).min() >= 0.0);
      CHECK(mRandomGenerator(i).max() <= 1.0);

      // Gather statistics
#pragma omp critical
      {
        minNumFlaws = min(minNumFlaws, size_t(mNumFlaws(i)));
        maxNumFlaws = max(maxNumFlaws, size_t(mNumFlaws(i)));
        totalNumFlaws += size_t(mNumFlaws(i));
        epsMin = std::min(epsMin, mMinFlaw(i));
        epsMax = std::max(epsMax, mMaxFlaw(i));
        sumFlaws += sumFlawsi;
      }
    }
  }

  // Some diagnostic output.
  if (nused_global > 0) {
    minNumFlaws = allReduce(minNumFlaws, MPI_MIN, Communicator::communicator());
    maxNumFlaws = allReduce(maxNumFlaws, MPI_MAX, Communicator::communicator());
    totalNumFlaws = allReduce(totalNumFlaws, MPI_SUM, Communicator::communicator());
    epsMin = allReduce(epsMin, MPI_MIN, Communicator::communicator());
    epsMax = allReduce(epsMax, MPI_MAX, Communicator::communicator());
    sumFlaws = allReduce(sumFlaws, MPI_SUM, Communicator::communicator());
    if (Process::getRank() == 0) {
      cerr << "ProbabilisticDamageModel for " << nodes.name() << ":" << endl
           << " Min, max, max/min volumes: " << mVmin << " " << mVmax << " " << mVmax*safeInv(mVmin) << endl
           << "    Min num flaws per node: " << minNumFlaws << endl
           << "    Max num flaws per node: " << maxNumFlaws << endl
           << "    Total num flaws       : " << totalNumFlaws << endl
           << "    Avg flaws per node    : " << totalNumFlaws / std::max(1, nused_global) << endl
           << "    Min flaw strain       : " << epsMin << endl
           << "    Max flaw strain       : " << epsMax << endl
           << "    Avg node failure      : " << sumFlaws / std::max(1, nused_global) << endl;
    }
  }
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
//
// In this model we compute the scalar damage derivative assuming unresolved
// crack growth for every point. However, that is not applied in the tensor
// damage update policy unless the flaws are actually activated.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
evaluateDerivatives(const Scalar time,
                    const Scalar dt,
                    const DataBase<Dimension>& dataBase,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Set the scalar magnitude of the damage evolution.
  const auto* nodeListPtr = &(this->nodeList());
  auto&       DDDt = derivs.field(state.buildFieldKey(ProbabilisticDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeListPtr->name()), 0.0);
  this->computeScalarDDDt(dataBase,
                          state,
                          time,
                          dt,
                          DDDt);
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename ProbabilisticDamageModel<Dimension>::TimeStepType
ProbabilisticDamageModel<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {

  // // Look at how quickly we're trying to change the damage.
  // double dt = DBL_MAX;
  // const Field<Dimension, SymTensor>& damage = this->nodeList().damage();
  // const ConnectivityMap<Dimension>& connectivityMap = dataBase.connectivityMap();
  // const vector<const NodeList<Dimension>*>& nodeLists = connectivityMap.nodeLists();
  // const size_t nodeListi = distance(nodeLists.begin(), find(nodeLists.begin(), nodeLists.end(), &(this->nodeList())));
  // for (typename ConnectivityMap<Dimension>::const_iterator iItr = connectivityMap.begin(nodeListi);
  //      iItr != connectivityMap.end(nodeListi);
  //      ++iItr) {
  //   const int i = *iItr;
  //   const double D0 = damage(i).Trace() / Dimension::nDim;
  //   dt = min(dt, 0.8*max(D0, 1.0 - D0)/
  //            std::sqrt(mDdamageDt(i)*mDdamageDt(i) + 1.0e-20));
  // }
  // return TimeStepType(dt, "Rate of damage change");

  return TimeStepType(1.0e100, "Rate of damage change -- NO VOTE.");
}

//------------------------------------------------------------------------------
// Register our state.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
registerState(DataBase<Dimension>& dataBase,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Register Youngs modulus and the longitudinal sound speed.
  PolicyPointer EPolicy(new YoungsModulusPolicy<Dimension>());
  PolicyPointer clPolicy(new LongitudinalSoundSpeedPolicy<Dimension>());
  state.enroll(mYoungsModulus, EPolicy);
  state.enroll(mLongitudinalSoundSpeed, clPolicy);

  // Set the initial values for the Youngs modulus, sound speed, and pressure.
  typename StateDerivatives<Dimension>::PackageList dummyPackages;
  StateDerivatives<Dimension> derivs(dataBase, dummyPackages);
  EPolicy->update(state.key(mYoungsModulus), state, derivs, 1.0, 0.0, 0.0);
  clPolicy->update(state.key(mLongitudinalSoundSpeed), state, derivs, 1.0, 0.0, 0.0);

  // Register the strain and effective strain.
  PolicyPointer effectiveStrainPolicy(new TensorStrainPolicy<Dimension>(mStrainAlgorithm));
  state.enroll(mStrain);
  state.enroll(mEffectiveStrain, effectiveStrainPolicy);

  // Register the damage and state it requires.
  // Note we are overriding the default no-op policy for the damage
  // as originally registered by the SolidSPHHydroBase class.
  auto& damage = this->nodeList().damage();
  PolicyPointer damagePolicy(new ProbabilisticDamagePolicy<Dimension>(mDamageInCompression,
                                                                      mkWeibull,
                                                                      mmWeibull,
                                                                      mMinFlawsPerNode,
                                                                      mVmin,
                                                                      mVmax));
  state.enroll(damage, damagePolicy);
  state.enroll(mNumFlaws);
  state.enroll(mMinFlaw);
  state.enroll(mMaxFlaw);
  state.enroll(mInitialVolume);
  state.enroll(mRandomGenerator);

  // Mask out nodes beyond the critical damage threshold from setting the timestep.
  auto maskKey = state.buildFieldKey(HydroFieldNames::timeStepMask, this->nodeList().name());
  auto& mask = state.field(maskKey, 0);
  const auto nlocal = this->nodeList().numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    if (damage(i).Trace() > mCriticalDamageThreshold) mask(i) = 0;
  }
}

//------------------------------------------------------------------------------
// Register the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mDdamageDt);
}

//------------------------------------------------------------------------------
// Apply the boundary conditions to the ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  CHECK(state.registered(DKey));
  auto& D = state.field(DKey, SymTensor::zero);

  // Apply ghost boundaries to the damage.
  for (auto boundaryItr = this->boundaryBegin();
       boundaryItr < this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyGhostBoundary(D);
  }
}

//------------------------------------------------------------------------------
// Enforce boundary conditions for the physics specific fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  CHECK(state.registered(DKey));
  auto& D = state.field(DKey, SymTensor::zero);

  // Enforce!
  for (auto boundaryItr = this->boundaryBegin(); 
       boundaryItr < this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceBoundary(D);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DamageModel<Dimension>::dumpState(file, pathName);
  file.write(mNumFlaws, pathName + "/numFlaws");
  file.write(mMinFlaw, pathName + "/minFlaw");
  file.write(mMaxFlaw, pathName + "/maxFlaw");
  file.write(mStrain, pathName + "/strain");
  file.write(mEffectiveStrain, pathName + "/effectiveStrain");
  file.write(mDdamageDt, pathName + "/DdamageDt");
  file.write(mMask, pathName + "/mask");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
restoreState(const FileIO& file, const string& pathName) {
  DamageModel<Dimension>::restoreState(file, pathName);
  file.read(mNumFlaws, pathName + "/numFlaws");
  file.read(mMinFlaw, pathName + "/minFlaw");
  file.read(mMaxFlaw, pathName + "/maxFlaw");
  file.read(mStrain, pathName + "/strain");
  file.read(mEffectiveStrain, pathName + "/effectiveStrain");
  file.read(mDdamageDt, pathName + "/DdamageDt");
  file.read(mMask, pathName + "/mask");
}

}

