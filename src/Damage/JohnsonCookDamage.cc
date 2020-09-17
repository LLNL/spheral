//---------------------------------Spheral++----------------------------------//
// JohnsonCookDamage -- an implementation of a Johnson-Cook damage law.
//
// Created by JMO, Mon Jul  9 08:21:23 PDT 2018
//----------------------------------------------------------------------------//
#include "FileIO/FileIO.hh"
#include "JohnsonCookDamage.hh"
#include "JohnsonCookFailureStrainPolicy.hh"
#include "JohnsonCookDamagePolicy.hh"
#include "EffectiveTensorDamagePolicy.hh"
#include "Strength/SolidFieldNames.hh"
#include "Strength/MeltEnergyPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Field/FieldList.hh"
#include "Boundary/Boundary.hh"
#include "Neighbor/Neighbor.hh"
#include "Utilities/globalNodeIDs.hh"

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
JohnsonCookDamage<Dimension>::
JohnsonCookDamage(SolidNodeList<Dimension>& nodeList,
                  const Field<Dimension, Scalar>& D1,
                  const Field<Dimension, Scalar>& D2,
                  const double D3,
                  const double D4,
                  const double D5,
                  const double epsilondot0,
                  const double Tcrit,
                  const double sigmamax,
                  const double efailmin):
  mNodeList(nodeList),
  mD1("D1_" + nodeList.name(), D1),
  mD2("D2_" + nodeList.name(), D2),
  mFailureStrain(SolidFieldNames::flaws, nodeList),
  mMeltSpecificEnergy(SolidFieldNames::meltSpecificEnergy, nodeList),
  mNewEffectiveDamage(EffectiveTensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::effectiveTensorDamage, nodeList),
  mD3(D3),
  mD4(D4),
  mD5(D5),
  mepsilondot0(epsilondot0),
  mTcrit(Tcrit),
  msigmamax(sigmamax),
  mefailmin(efailmin),
  mRestart(registerWithRestart(*this)) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookDamage<Dimension>::
~JohnsonCookDamage() {
}

//------------------------------------------------------------------------------
// Increment the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamage<Dimension>::
evaluateDerivatives(const Scalar /*time*/,
                    const Scalar /*dt*/,
                    const DataBase<Dimension>& /*dataBase*/,
                    const State<Dimension>& state,
                    StateDerivatives<Dimension>& derivs) const {

  // Fow now we update the effective damage as simply a copy of the damage.
  const auto& D = state.field(state.buildFieldKey(SolidFieldNames::tensorDamage, mNodeList.name()), SymTensor::zero);
  auto& Deff = derivs.field(state.buildFieldKey(EffectiveTensorDamagePolicy<Dimension>::prefix() + SolidFieldNames::effectiveTensorDamage, mNodeList.name()), SymTensor::zero);
  Deff = D;
}

//------------------------------------------------------------------------------
// Vote on a time step.
//------------------------------------------------------------------------------
template<typename Dimension>
typename JohnsonCookDamage<Dimension>::TimeStepType
JohnsonCookDamage<Dimension>::
dt(const DataBase<Dimension>& /*dataBase*/, 
   const State<Dimension>& /*state*/,
   const StateDerivatives<Dimension>& /*derivs*/,
   const Scalar /*currentTime*/) const {
  return TimeStepType(1.0e100, "Rate of damage change -- NO VOTE.");
}

//------------------------------------------------------------------------------
// Register the state and update policies.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamage<Dimension>::
registerState(DataBase<Dimension>& /*dataBase*/,
              State<Dimension>& state) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Register the failure strains for updating.
  PolicyPointer flawPolicy(new JohnsonCookFailureStrainPolicy<Dimension>(mD1,
                                                                         mD2,
                                                                         mD3,
                                                                         mD4,
                                                                         mD5,
                                                                         mepsilondot0,
                                                                         msigmamax,
                                                                         mefailmin,
                                                                         mTcrit));
  state.enroll(mFailureStrain, flawPolicy);

  // Register the damage and effective damage.
  PolicyPointer damagePolicy(new JohnsonCookDamagePolicy<Dimension>());
  PolicyPointer effDamagePolicy(new EffectiveTensorDamagePolicy<Dimension>());
  state.enroll(mNodeList.damage(), damagePolicy);
  state.enroll(mNodeList.effectiveDamage(), effDamagePolicy);

  // We also require the melt energy.
  PolicyPointer meltPolicy(new MeltEnergyPolicy<Dimension>());
  state.enroll(mMeltSpecificEnergy, meltPolicy);
}

//------------------------------------------------------------------------------
// Register the derivatives.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamage<Dimension>::
registerDerivatives(DataBase<Dimension>& /*dataBase*/,
                    StateDerivatives<Dimension>& derivs) {
  derivs.enroll(mNewEffectiveDamage);
}

//------------------------------------------------------------------------------
// Apply the boundary conditions to the ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamage<Dimension>::
applyGhostBoundaries(State<Dimension>& state,
                     StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  const Key DeffKey = state.buildFieldKey(SolidFieldNames::effectiveTensorDamage, nodeListName);
  CHECK(state.registered(DKey));
  CHECK(state.registered(DeffKey));
  auto& D = state.field(DKey, SymTensor::zero);
  auto& Deff = state.field(DeffKey, SymTensor::zero);

  // Apply ghost boundaries to the damage.
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin();
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->applyGhostBoundary(D);
    (*boundaryItr)->applyGhostBoundary(Deff);
  }
}

//------------------------------------------------------------------------------
// Enforce boundary conditions for the physics specific fields.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamage<Dimension>::
enforceBoundaries(State<Dimension>& state,
                  StateDerivatives<Dimension>& /*derivs*/) {

  // Grab this models damage field from the state.
  typedef typename State<Dimension>::KeyType Key;
  const Key nodeListName = this->nodeList().name();
  const Key DKey = state.buildFieldKey(SolidFieldNames::tensorDamage, nodeListName);
  const Key DeffKey = state.buildFieldKey(SolidFieldNames::effectiveTensorDamage, nodeListName);
  const Key gradDKey = state.buildFieldKey(SolidFieldNames::damageGradient, nodeListName);
  CHECK(state.registered(DKey));
  CHECK(state.registered(DeffKey));
  CHECK(state.registered(gradDKey));
  auto& D = state.field(DKey, SymTensor::zero);
  auto& Deff = state.field(DeffKey, SymTensor::zero);

  // Enforce!
  for (ConstBoundaryIterator boundaryItr = this->boundaryBegin(); 
       boundaryItr != this->boundaryEnd();
       ++boundaryItr) {
    (*boundaryItr)->enforceBoundary(D);
    (*boundaryItr)->enforceBoundary(Deff);
  }
}

//------------------------------------------------------------------------------
// Dump the current state to the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamage<Dimension>::
dumpState(FileIO& file, const string& pathName0) const {
  const string pathName = pathName0 + "/" + mNodeList.name();
  file.write(mD1, pathName + "/D1");
  file.write(mD2, pathName + "/D2");
  file.write(mFailureStrain, pathName + "/failureStrain");
  file.write(mMeltSpecificEnergy, pathName + "/meltSpecificEnergy");
}

//------------------------------------------------------------------------------
// Restore the state from the given file.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookDamage<Dimension>::
restoreState(const FileIO& file, const string& pathName0) {
  const string pathName = pathName0 + "/" + mNodeList.name();
  file.read(mD1, pathName + "/D1");
  file.read(mD2, pathName + "/D2");
  file.read(mFailureStrain, pathName + "/failureStrain");
  file.read(mMeltSpecificEnergy, pathName + "/meltSpecificEnergy");
}

}
