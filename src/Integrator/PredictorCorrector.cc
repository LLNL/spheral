//---------------------------------Spheral++----------------------------------//
// Second order predictor corrector time integrator.
// Based on the predictor corrector scheme described in Monaghans tensile
// instability paper (Monaghan 2000, JCP 159, 290.)
//
// Created by JMO, Wed Dec  4 21:39:36 PST 2002
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "PredictorCorrector.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"

#include "DBC.hh"

namespace Spheral {
namespace IntegratorSpace {

using namespace std;

using DataBaseSpace::DataBase;
using FieldSpace::FieldList;
using PhysicsSpace::Physics;

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>::PredictorCorrector():
  Integrator<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>::
PredictorCorrector(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>::
PredictorCorrector(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>::~PredictorCorrector() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>&
PredictorCorrector<Dimension>::
operator=(const PredictorCorrector<Dimension>& rhs) {
  if (this != &rhs) {
    Integrator<Dimension>::operator=(rhs);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PredictorCorrector<Dimension>::
step(typename Dimension::Scalar maxTime) {

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Construct the state and derivatives, and initalize the integrator.
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  this->initialize(state, derivs);

  // Determine the minimum timestep across all packages.
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);

  // Zero out the derivatives.
  derivs.Zero();

  // Evaluate the beginning of step derivatives.
  this->evaluateDerivatives(t, dt, db, state, derivs);
  this->finalizeDerivatives(t, dt, db, state, derivs);

  // Copy the beginning of step state and derivatives.
  State<Dimension> state0(state);
  StateDerivatives<Dimension> derivs0(derivs);
  state0.copyState();
  derivs0.copyState();

  // Predictor step -- trial advance to the end of step.
  state.update(derivs, dt, t, dt);

  // Enforce Boundary conditions.
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
                               
  // Do any physics specific stuff relating to the fact the state was just updated.
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Loop over the physics packages and perform any necessary initializations.
  this->preStepInitialize(t + dt, dt, state, derivs);

  // Zero out the stored derivatives.
  derivs.Zero();

  // Now loop over the packages and evaluate the derivatives at the
  // predicted end of step state.
  this->evaluateDerivatives(t + dt, dt, db, state, derivs);
  this->finalizeDerivatives(t + dt, dt, db, state, derivs);

  // Apply the corrector stage.
  const Scalar hdt = 0.5*dt;
  this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs0, hdt, t, hdt);
  this->applyGhostBoundaries(state, derivs);
  state.update(derivs, hdt, t + hdt, hdt);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);

  // Do any physics specific stuff relating to the fact the state was just updated.
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Apply any physics specific finalizations.
  this->finalize(t + dt, dt, state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
}
}
}

