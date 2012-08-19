//---------------------------------Spheral++----------------------------------//
// CheapSynchronousRK2 -- Advance the set of Physics packages in time using 
// my standard cheat on the second order Runge-Kutta.  We just skip the end of
// step derivative evaluation, and assume that we can use the previous cycles
// mid-step derivatives to determine this steps mid-step state.
//
// Created by JMO, Mon Jun 12 17:52:27 PDT 2000
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "CheapSynchronousRK2.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"

#include "DBC.hh"
#include "cdebug.hh"

#include "TAU.h"

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
CheapSynchronousRK2<Dimension>::CheapSynchronousRK2():
  Integrator<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
CheapSynchronousRK2<Dimension>::
CheapSynchronousRK2(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
CheapSynchronousRK2<Dimension>::
CheapSynchronousRK2(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
CheapSynchronousRK2<Dimension>::~CheapSynchronousRK2() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
CheapSynchronousRK2<Dimension>&
CheapSynchronousRK2<Dimension>::
operator=(const CheapSynchronousRK2<Dimension>& rhs) {
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
CheapSynchronousRK2<Dimension>::
step(typename Dimension::Scalar maxTime) {

  // TAU timers.
  TAU_PROFILE("CheapSynchronousRK2", "::step", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2ConstructState, "CheapSynchronousRK2", "::step : Construct state ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Dt, "CheapSynchronousRK2", "::step : Set dt", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2CopyFields, "CheapSynchronousRK2", "::step : Copy initial state", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2MidStep, "CheapSynchronousRK2", "::step : Advance to mid-step", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Boundaries1, "CheapSynchronousRK2", "::step: Set boundaries 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2PostState1, "CheapSynchronousRK2", "::step: Post state update 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2StepInitialize1, "CheapSynchronousRK2", "::step : Pre-step initialize 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Zero2, "CheapSynchronousRK2", "::step : Zero derivs  2          ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2EvalDerivs2, "CheapSynchronousRK2", "::step : Eval derivs mid step    ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2EndStep, "CheapSynchronousRK2", "::step : Advance to end of step  ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Boundaries2, "CheapSynchronousRK2", "::step : Set boundaries 2        ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2PostState2, "CheapSynchronousRK2", "::step: Post state update 2", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Finalize, "CheapSynchronousRK2", "::step : finalize physics", TAU_USER);

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Construct the state and derivatives, and initalize the integrator.
  TAU_PROFILE_START(TimeRK2ConstructState);
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  state.timeAdvanceOnly(true);
  this->initialize(state, derivs);
  TAU_PROFILE_STOP(TimeRK2ConstructState);

  // Determine the minimum timestep across all packages.
  TAU_PROFILE_START(TimeRK2Dt);
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);
  TAU_PROFILE_STOP(TimeRK2Dt);

  // Copy the beginning of step state.
  TAU_PROFILE_START(TimeRK2CopyFields);
  State<Dimension> state0(state);
  state0.copyState();
  state0.timeAdvanceOnly(false);
  TAU_PROFILE_STOP(TimeRK2CopyFields);

  // Trial advance the state to the mid timestep point.
  TAU_PROFILE_START(TimeRK2MidStep);
  const double hdt = 0.5*dt;
  state.update(derivs, hdt, t, hdt);
  this->currentTime(t + hdt);
  TAU_PROFILE_STOP(TimeRK2MidStep);

  // Enforce Boundary conditions.
  TAU_PROFILE_START(TimeRK2Boundaries1);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  TAU_PROFILE_STOP(TimeRK2Boundaries1);
                                  
  // Do any physics specific stuff relating to the fact the state was just updated.
  TAU_PROFILE_START(TimeRK2PostState1);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();
  TAU_PROFILE_STOP(TimeRK2PostState1);

  // Loop over the physics packages and perform any necessary initializations.
  TAU_PROFILE_START(TimeRK2StepInitialize1);
  this->preStepInitialize(t + hdt, hdt, state, derivs);
  TAU_PROFILE_STOP(TimeRK2StepInitialize1);

  // Zero out the stored derivatives.
  TAU_PROFILE_START(TimeRK2Zero2);
  derivs.Zero();
  TAU_PROFILE_STOP(TimeRK2Zero2);

  // Now loop over the packages and evaluate the derivatives at the
  // midpoint.
  TAU_PROFILE_START(TimeRK2EvalDerivs2);
  this->evaluateDerivatives(t + hdt, hdt, db, state, derivs);
  this->finalizeDerivatives(t + hdt, hdt, db, state, derivs);
  TAU_PROFILE_STOP(TimeRK2EvalDerivs2);

  // Advance the state from the beginning of the cycle using the midpoint 
  // derivatives.
  TAU_PROFILE_START(TimeRK2EndStep);
  state.timeAdvanceOnly(false);
  this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);
  this->currentTime(t + dt);
  TAU_PROFILE_STOP(TimeRK2EndStep);

  // Enforce boundaries.
  TAU_PROFILE_START(TimeRK2Boundaries2);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  TAU_PROFILE_STOP(TimeRK2Boundaries2);

  // Do any physics specific stuff relating to the fact the state was just updated.
  TAU_PROFILE_START(TimeRK2PostState2);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();
  TAU_PROFILE_STOP(TimeRK2PostState2);

  // Apply any physics specific finalizations.
  TAU_PROFILE_START(TimeRK2Finalize);
  this->finalize(t + dt, dt, state, derivs);
  TAU_PROFILE_STOP(TimeRK2Finalize);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
}
}
}

