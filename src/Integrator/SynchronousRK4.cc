//---------------------------------Spheral++----------------------------------//
// SynchronousRK4 -- Advance the set of Physics packages in time using fourth
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun 14 10:22:47 PDT 2012
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "SynchronousRK4.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"
#include "DBC.hh"

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
SynchronousRK4<Dimension>::SynchronousRK4():
  Integrator<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK4<Dimension>::
SynchronousRK4(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK4<Dimension>::
SynchronousRK4(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK4<Dimension>::~SynchronousRK4() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK4<Dimension>&
SynchronousRK4<Dimension>::
operator=(const SynchronousRK4<Dimension>& rhs) {
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
SynchronousRK4<Dimension>::
step(typename Dimension::Scalar maxTime) {

  // TAU timers.
  TAU_PROFILE("SynchronousRK4", "::step", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4ConstructState, "SynchronousRK4", "::step : Construct state ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4Dt, "SynchronousRK4", "::step : Set dt", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4Zero, "SynchronousRK4", "::step : Zero derivs             ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4EvalDerivs1, "SynchronousRK4", "::step : Eval derivs @ t0  ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4CopyState, "SynchronousRK4", "::step : Copy initial state", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4Stage1, "SynchronousRK4", "::step : stage 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4Stage2, "SynchronousRK4", "::step : stage 2", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4Stage3, "SynchronousRK4", "::step : stage 3", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4Stage4, "SynchronousRK4", "::step : stage 4", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4AdvanceState, "SynchronousRK4", "::step : Advance to end of step  ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK4Finalize, "SynchronousRK4", "::step : finalize physics", TAU_USER);

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Construct the state and derivatives, and initalize the integrator.
  TAU_PROFILE_START(TimeRK4ConstructState);
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs1(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  this->initialize(state, derivs1);
  TAU_PROFILE_STOP(TimeRK4ConstructState);

  // Determine the minimum timestep across all packages.
  TAU_PROFILE_START(TimeRK4Dt);
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs1);
  TAU_PROFILE_STOP(TimeRK4Dt);

  // Zero out the derivatives, and make some independent copies
  TAU_PROFILE_START(TimeRK4Zero);
  derivs1.Zero();
  StateDerivatives<Dimension> derivs2(derivs1), derivs3(derivs1), derivs4(derivs1);
  derivs2.copyState();
  derivs3.copyState();
  derivs4.copyState();
  TAU_PROFILE_STOP(TimeRK4Zero);

  // Make a copy of the state we'll use for our intermediate estimates.
  TAU_PROFILE_START(TimeRK4CopyState);
  State<Dimension> tmpstate(state);
  tmpstate.copyState();
  TAU_PROFILE_STOP(TimeRK4CopyState);

  // Stage 1:
  // Get derivs1(t_n, state(t_n))
  TAU_PROFILE_START(TimeRK4Stage1);
  this->evaluateDerivatives(t, dt, db, state, derivs1);
  this->finalizeDerivatives(t, dt, db, state, derivs1);
  TAU_PROFILE_STOP(TimeRK4Stage1);

  // Stage 2: 
  // Get derivs2(t_n + 0.5*dt, state(t_n + 0.5*dt*derivs1))
  TAU_PROFILE_START(TimeRK4Stage2);
  tmpstate.update(derivs1, 0.5*dt, t, 0.5*dt);
  this->enforceBoundaries(tmpstate, derivs1);
  this->applyGhostBoundaries(tmpstate, derivs1);
  this->postStateUpdate(db, tmpstate, derivs1);
  this->finalizeGhostBoundaries();
  this->preStepInitialize(t + 0.5*dt, 0.5*dt, tmpstate, derivs2);
  this->evaluateDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs2);
  this->finalizeDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs2);
  TAU_PROFILE_STOP(TimeRK4Stage2);

  // Stage 3: 
  // Get derivs3(t_n + 0.5*dt, state(t_n + 0.5*dt*derivs2))
  TAU_PROFILE_START(TimeRK4Stage3);
  tmpstate = state;
  tmpstate.copyState();
  tmpstate.update(derivs2, 0.5*dt, t, 0.5*dt);
  this->enforceBoundaries(tmpstate, derivs2);
  this->applyGhostBoundaries(tmpstate, derivs2);
  this->postStateUpdate(db, tmpstate, derivs2);
  this->finalizeGhostBoundaries();
  this->preStepInitialize(t + 0.5*dt, 0.5*dt, tmpstate, derivs3);
  this->evaluateDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs3);
  this->finalizeDerivatives(t + 0.5*dt, 0.5*dt, db, tmpstate, derivs3);
  TAU_PROFILE_STOP(TimeRK4Stage3);

  // Stage 4: 
  // Get derivs3(t_n + dt, state(t_n + dt*derivs3))
  TAU_PROFILE_START(TimeRK4Stage4);
  tmpstate = state;
  tmpstate.copyState();
  tmpstate.update(derivs3, dt, t, dt);
  this->enforceBoundaries(tmpstate, derivs3);
  this->applyGhostBoundaries(tmpstate, derivs3);
  this->postStateUpdate(db, tmpstate, derivs3);
  this->finalizeGhostBoundaries();
  this->preStepInitialize(t + dt, dt, tmpstate, derivs4);
  this->evaluateDerivatives(t + dt, dt, db, tmpstate, derivs4);
  this->finalizeDerivatives(t + dt, dt, db, tmpstate, derivs4);
  TAU_PROFILE_STOP(TimeRK4Stage4);

  // Advance.
  // Now we can apply the RK4 algorithm to advance the actual state over the full time step.
  // Conceptually we are doing:
  //   state(t_n + dt) = state(t_n) + dt/6*(derivs1 + 2*derivs2 + 2*derivs3 + derivs4)
  TAU_PROFILE_START(TimeRK4AdvanceState);
  state.update(derivs1, dt/6.0, t, dt);
  state.update(derivs2, dt/3.0, t, dt);
  state.update(derivs3, dt/3.0, t, dt);
  state.update(derivs4, dt/6.0, t, dt);
  this->enforceBoundaries(state, derivs4);
  this->applyGhostBoundaries(state, derivs4);
  this->postStateUpdate(db, state, derivs4);
  this->finalizeGhostBoundaries();
  TAU_PROFILE_STOP(TimeRK4AdvanceState);

  // Apply any physics specific finalizations.
  TAU_PROFILE_START(TimeRK4Finalize);
  this->finalize(t + dt, dt, state, derivs4);
  TAU_PROFILE_STOP(TimeRK4Finalize);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
}
}
}

