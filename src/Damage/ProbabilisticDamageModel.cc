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
                         const DamageCouplingAlgorithm damageCouplingAlgorithm,
                         const TensorStrainAlgorithm strainAlgorithm,
                         const bool damageInCompression):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, damageCouplingAlgorithm),
  mStrainAlgorithm(strainAlgorithm),
  mDamageInCompression(damageInCompression),
  mkWeibull(kWeibull),
  mmWeibull(mWeibull),
  mVmin(0.0),
  mVmax(0.0),
  mSeed(seed),
  mMinFlawsPerNode(minFlawsPerNode),
  mNumFlaws(SolidFieldNames::numFlaws, nodeList),
  mNumFlawsActivated(SolidFieldNames::numFlawsActivated, nodeList),
  mCurrentFlaw(SolidFieldNames::currentFlaw, nodeList, std::numeric_limits<double>::max()),
  mInitialVolume(SolidFieldNames::initialVolume, nodeList),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mDdamageDt(ProbabilisticDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList),
  mRandomGenerators(SolidFieldNames::randomGenerators, nodeList) {
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
  const auto nlocal = nodes.numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    CHECK(mass(i) > 0.0 and rho(i) > 0.0);
    mInitialVolume(i) = mass(i)/rho(i);
    Key seedi = mSeed;
    boost::hash_combine(seedi, keys(i));
    mRandomGenerators(i).seed(seedi);
    mRandomGenerators(i)();                // Recommended to discard first value in sequence
  }

  // Find the minimum and maximum node volumes.
  mVmin = mInitialVolume.min();
  mVmax = mInitialVolume.max();
  CHECK(mVmin > 0.0);
  CHECK(mVmax >= mVmin);

  // Compute the maximum strain we expect for the minimum volume.
  const auto epsMax2m = mMinFlawsPerNode/(mkWeibull*mVmin);  // epsmax ** m

  // Based on this compute the maximum number of flaws any node will have.  We'll use this to
  // spin the random number generator without extra communiction.
  const auto maxFlawsPerNode = mMinFlawsPerNode*std::max(1u, unsigned(mkWeibull*mVmax*epsMax2m + 0.5));

  // Now find the minimum actiavation strain flaw each point will start with.
  const auto mInv = 1.0/mmWeibull;
#pragma omp parallel for
  for (auto i = 0u; i < nlocal; ++i) {
    mNumFlaws(i) = std::max(size_t(1), std::min(maxFlawsPerNode, size_t(mkWeibull*mInitialVolume(i)*epsMax2m + 0.5)));
    const auto Ai = mNumFlaws(i)/(mkWeibull*mInitialVolume(i));
    CHECK(Ai > 0.0);
    for (auto j = 0u; j < mNumFlaws(i); ++j) {
      mCurrentFlaw(i) = std::min(mCurrentFlaw(i), pow(Ai * mRandomGenerators(i)(), mInv));
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
  PolicyPointer damagePolicy(new ProbabilisticDamagePolicy<Dimension>(mDamageInCompression,
                                                                      mkWeibull,
                                                                      mmWeibull,
                                                                      mMinFlawsPerNode,
                                                                      mVmin,
                                                                      mVmax));
  state.enroll(this->nodeList().damage(), damagePolicy);
  state.enroll(mNumFlaws);
  state.enroll(mNumFlawsActivated);
  state.enroll(mCurrentFlaw);
  state.enroll(mInitialVolume);
  state.enroll(mRandomGenerators);
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
  file.write(mNumFlawsActivated, pathName + "/numFlawsActivated");
  file.write(mCurrentFlaw, pathName + "/currentFlaw");
  file.write(mStrain, pathName + "/strain");
  file.write(mEffectiveStrain, pathName + "/effectiveStrain");
  file.write(mDdamageDt, pathName + "/DdamageDt");
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
  file.read(mNumFlawsActivated, pathName + "/numFlawsActivated");
  file.read(mCurrentFlaw, pathName + "/currentFlaw");
  file.read(mStrain, pathName + "/strain");
  file.read(mEffectiveStrain, pathName + "/effectiveStrain");
  file.read(mDdamageDt, pathName + "/DdamageDt");
}

}

