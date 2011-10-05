//---------------------------------Spheral++----------------------------------//
// SynchronousRK2 -- Advance the set of Physics packages in time using second
// order Runge-Kutta.  All packages are advanced at one timestep simultaneously
// each step, i.e., synchronously.
//
// Created by JMO, Thu Jun  1 22:26:45 PDT 2000
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "SynchronousRK2.hh"
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
SynchronousRK2<Dimension>::SynchronousRK2():
  Integrator<Dimension>() {
  cdebug << "SynchronousRK2::SynchronousRK2()" << endl;
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK2<Dimension>::
SynchronousRK2(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
  cdebug << "SynchronousRK2::SynchronousRK2(DataBase)" << endl;
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK2<Dimension>::
SynchronousRK2(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
  cdebug << "SynchronousRK2::SynchronousRK2(DataBase, PhysicsPackages)" << endl;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK2<Dimension>::~SynchronousRK2() {
  cdebug << "SynchronousRK2::~SynchronousRK2()" << endl;
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK2<Dimension>&
SynchronousRK2<Dimension>::
operator=(const SynchronousRK2<Dimension>& rhs) {
  cdebug << "SynchronousRK2::operator=()" << endl;
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
SynchronousRK2<Dimension>::
step(typename Dimension::Scalar maxTime) {

  // TAU timers.
  TAU_PROFILE("SynchronousRK2", "::step", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2ConstructState, "SynchronousRK2", "::step : Construct state ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Dt, "SynchronousRK2", "::step : Set dt", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Zero1, "SynchronousRK2", "::step : Zero derivs  1          ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2EvalDerivs1, "SynchronousRK2", "::step : Eval derivs @ t0  ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2CopyFields, "SynchronousRK2", "::step : Copy initial state", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2MidStep, "SynchronousRK2", "::step : Advance to mid-step", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Boundaries1, "SynchronousRK2", "::step: Set boundaries 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2PostState1, "SynchronousRK2", "::step: Post state update 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2StepInitialize1, "SynchronousRK2", "::step : Pre-step initialize 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Zero2, "SynchronousRK2", "::step : Zero derivs  2          ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2EvalDerivs2, "SynchronousRK2", "::step : Eval derivs mid step    ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2EndStep, "SynchronousRK2", "::step : Advance to end of step  ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Boundaries2, "SynchronousRK2", "::step : Set boundaries 2        ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2PostState2, "SynchronousRK2", "::step: Post state update 2", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK2Finalize, "SynchronousRK2", "::step : finalize physics", TAU_USER);

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Construct the state and derivatives, and initalize the integrator.
  TAU_PROFILE_START(TimeRK2ConstructState);
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  initialize(state, derivs);
  TAU_PROFILE_STOP(TimeRK2ConstructState);

  // Determine the minimum timestep across all packages.
  TAU_PROFILE_START(TimeRK2Dt);
  const Scalar dt = selectDt(min(this->dtMin(), maxTime - t),
                             min(this->dtMax(), maxTime - t),
                             state,
                             derivs);
  cdebug << "SynchronousRK2::step: chose dt = " << dt << endl;
  TAU_PROFILE_STOP(TimeRK2Dt);

  // Zero out the derivatives.
  TAU_PROFILE_START(TimeRK2Zero1);
  derivs.Zero();
  TAU_PROFILE_STOP(TimeRK2Zero1);

  // Evaluate the beginning of step derivatives.
  TAU_PROFILE_START(TimeRK2EvalDerivs1);
  const double hdt = 0.5*dt;
  this->evaluateDerivatives(t, hdt, db, state, derivs);
  this->finalizeDerivatives(t, hdt, db, state, derivs);
  TAU_PROFILE_STOP(TimeRK2EvalDerivs1);

  // Copy the beginning of step state.
  TAU_PROFILE_START(TimeRK2CopyFields);
  State<Dimension> state0(state);
  state0.copyState();
  TAU_PROFILE_STOP(TimeRK2CopyFields);

  // Trial advance the state to the mid timestep point.
  TAU_PROFILE_START(TimeRK2MidStep);
  state.update(derivs, hdt, t, hdt);
  TAU_PROFILE_STOP(TimeRK2MidStep);

  // Enforce Boundary conditions (as a side effect this updates the
  // neighboring).
  TAU_PROFILE_START(TimeRK2Boundaries1);
  enforceBoundaries(state, derivs);
  applyGhostBoundaries(state, derivs);
  TAU_PROFILE_STOP(TimeRK2Boundaries1);
                                  
  // Do any physics specific stuff relating to the fact the state was just updated.
  TAU_PROFILE_START(TimeRK2PostState1);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();
  TAU_PROFILE_STOP(TimeRK2PostState1);

  // Loop over the physics packages and perform any necessary initializations.
  TAU_PROFILE_START(TimeRK2StepInitialize1);
  preStepInitialize(t + hdt, hdt, state, derivs);
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
  copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);
  TAU_PROFILE_STOP(TimeRK2EndStep);

  // Enforce boundaries.
  TAU_PROFILE_START(TimeRK2Boundaries2);
  enforceBoundaries(state, derivs);
  applyGhostBoundaries(state, derivs);
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
  currentCycle(this->currentCycle() + 1);
  currentTime(t + dt);
  lastDt(dt);
}
}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
namespace IntegratorSpace {
template class SynchronousRK2< Dim<1> >;
template class SynchronousRK2< Dim<2> >;
template class SynchronousRK2< Dim<3> >;
}
}
