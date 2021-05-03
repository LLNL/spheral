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

#include <string>
#include <vector>
#include <algorithm>
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
                         const double crackGrowthMultiplier,
                         const DamageCouplingAlgorithm damageCouplingAlgorithm,
                         const TensorStrainAlgorithm strainAlgorithm,
                         const bool damageInCompression):
  DamageModel<Dimension>(nodeList, W, crackGrowthMultiplier, damageCouplingAlgorithm),
  mStrainAlgorithm(strainAlgorithm),
  mDamageInCompression(damageInCompression),
  mNumFlaws(SolidFieldNames::numFlaws, nodeList),
  mNumFlawsActivated(SolidFieldNames::numFlawsActivated, nodeList),
  mCurrentFlaw(SolidFieldNames::currentFlaw, nodeList),
  mYoungsModulus(SolidFieldNames::YoungsModulus, nodeList),
  mLongitudinalSoundSpeed(SolidFieldNames::longitudinalSoundSpeed, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mDdamageDt(ProbabilisticDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrainAlgorithm(strainAlgorithm) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ProbabilisticDamageModel<Dimension>::
~ProbabilisticDamageModel() {
}

//------------------------------------------------------------------------------
// Evaluate derivatives.
// For this model we evaluate the derivative of the tensor damage field.
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

  typedef typename State<Dimension>::KeyType Key;
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

  // Register the damage.
  // Note we are overriding the default no-op policy for the damage
  // as originally registered by the SolidSPHHydroBase class.
  PolicyPointer damagePolicy(new ProbabilisticDamagePolicy<Dimension>(*this));
  state.enroll(this->nodeList().damage(), damagePolicy);

  // Mask out nodes beyond the critical damage threshold from setting the timestep.
  Key maskKey = state.buildFieldKey(HydroFieldNames::timeStepMask, this->nodeList().name());
  Field<Dimension, int>& mask = state.field(maskKey, 0);
  const Field<Dimension, SymTensor>& damage = this->nodeList().damage();
  for (auto i = 0u; i < this->nodeList().numInternalNodes(); ++i) {
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
  Field<Dimension, SymTensor>& D = state.field(DKey, SymTensor::zero);

  // Apply ghost boundaries to the damage.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
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
  Field<Dimension, SymTensor>& D = state.field(DKey, SymTensor::zero);

  // Enforce!
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceBoundary(D);
  }
}

//------------------------------------------------------------------------------
// Cull the flaws on each node to the single weakest one.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
cullToWeakestFlaws() {
  const auto n = mFlaws.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    auto& flaws = mFlaws[i];
    if (flaws.size() > 0) {
      const auto maxVal = *max_element(flaws.begin(), flaws.end());
      flaws = vector<double>(maxVal);
    }
  }
}

//------------------------------------------------------------------------------
// Compute a Field with the sum of the activation energies per node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
ProbabilisticDamageModel<Dimension>::
sumActivationEnergiesPerNode() const {
  auto& nodeList = this->nodeList();
  const auto n = mFlaws.numInternalElements();
  Field<Dimension, Scalar> result("Sum activation energies", nodeList);
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto& flaws = mFlaws(i);
    for (const auto fij: flaws) {
      result(i) += fij;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Compute a Field with the number of flaws per node.
//------------------------------------------------------------------------------
template<typename Dimension>
Field<Dimension, typename Dimension::Scalar>
ProbabilisticDamageModel<Dimension>::
numFlawsPerNode() const {
  auto& nodeList = this->nodeList();
  const auto n = mFlaws.numInternalElements();
  Field<Dimension, Scalar> result("num flaws", nodeList);
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    result(i) = flawsForNode(i).size();
  }
  return result;
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ProbabilisticDamageModel<Dimension>::
dumpState(FileIO& file, const string& pathName) const {
  DamageModel<Dimension>::dumpState(file, pathName);
  file.write(mFlaws, pathName + "/flaws");
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
  file.read(mFlaws, pathName + "/flaws");
  file.read(mStrain, pathName + "/strain");
  file.read(mEffectiveStrain, pathName + "/effectiveStrain");
  file.read(mDdamageDt, pathName + "/DdamageDt");
}

}

