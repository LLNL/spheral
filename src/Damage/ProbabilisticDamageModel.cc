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
  mDdamageDt(ProbabilisticDamagePolicy<Dimension>::prefix() + SolidFieldNames::scalarDamage, nodeList),
  mStrain(SolidFieldNames::strainTensor, nodeList),
  mEffectiveStrain(SolidFieldNames::effectiveStrainTensor, nodeList) {
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

  // Register the strain, effective strain, and stuff needed by the damage update policy.
  PolicyPointer effectiveStrainPolicy(new TensorStrainPolicy<Dimension>(mStrainAlgorithm));
  state.enroll(mStrain);
  state.enroll(mEffectiveStrain, effectiveStrainPolicy);

  // Register the damage.
  // Note we are overriding the default no-op policy for the damage
  // as originally registered by the SolidSPHHydroBase class.
  PolicyPointer damagePolicy(new ProbabilisticDamagePolicy<Dimension>(mDamageInCompression));
  state.enroll(this->nodeList().damage(), damagePolicy);
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

