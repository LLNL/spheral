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
#include "Utilities/DBC.hh"
#include "Utilities/Timer.hh"
#include "caliper/cali.h"

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

// Declare the timers
extern Timer TIME_CheapRK2;
extern Timer TIME_CheapRK2PreInit;
extern Timer TIME_CheapRK2Dt;
extern Timer TIME_CheapRK2CopyState;
extern Timer TIME_CheapRK2MidStep;
extern Timer TIME_CheapRK2EvalDerivs;
extern Timer TIME_CheapRK2EndStep;
extern Timer TIME_CheapRK2Finalize;
extern Timer TIME_CheapRK2EnforceBound;

namespace Spheral {

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
bool
CheapSynchronousRK2<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  TIME_CheapRK2.start();
  CALI_MARK_BEGIN("TIME_CheapRK2");

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Initalize the integrator.
  TIME_CheapRK2PreInit.start();
  CALI_MARK_BEGIN("TIME_CheapRK2PreInit");
  this->preStepInitialize(state, derivs);
  this->initializeDerivatives(t, 0.0, state, derivs);
  TIME_CheapRK2PreInit.stop();
  CALI_MARK_END("TIME_CheapRK2PreInit");

  // Determine the minimum timestep across all packages.
  TIME_CheapRK2Dt.start();
  CALI_MARK_BEGIN("TIME_CheapRK2Dt");
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);
  const double hdt = 0.5*dt;
  TIME_CheapRK2Dt.stop();
  CALI_MARK_END("TIME_CheapRK2Dt");

  // Copy the beginning of step state.
  TIME_CheapRK2CopyState.start();
  CALI_MARK_BEGIN("TIME_CheapRK2CopyState");
  State<Dimension> state0(state);
  state0.copyState();
  TIME_CheapRK2CopyState.stop();
  CALI_MARK_END("TIME_CheapRK2CopyState");

  // Trial advance the state to the mid timestep point.
  TIME_CheapRK2MidStep.start();
  CALI_MARK_BEGIN("TIME_CheapRK2MidStep");
  state.timeAdvanceOnly(true);
  state.update(derivs, hdt, t, hdt);
  this->currentTime(t + hdt);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(t + hdt, hdt, db, state, derivs);
  this->finalizeGhostBoundaries();
  TIME_CheapRK2MidStep.stop();
  CALI_MARK_END("TIME_CheapRK2MidStep");

  // Evaluate the derivatives at the midpoint.
  TIME_CheapRK2EvalDerivs.start();
  CALI_MARK_BEGIN("TIME_CheapRK2EvalDerivs");
  this->initializeDerivatives(t + hdt, hdt, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t + hdt, hdt, db, state, derivs);
  this->finalizeDerivatives(t + hdt, hdt, db, state, derivs);
  TIME_CheapRK2EvalDerivs.stop();
  CALI_MARK_END("TIME_CheapRK2EvalDerivs");

  // Check if the timestep is still a good idea...
  if (this->allowDtCheck()) {
    const auto dtnew = this->selectDt(min(this->dtMin(), maxTime - t),
                                      min(this->dtMax(), maxTime - t),
                                      state,
                                      derivs);
    if (dtnew < 0.5*dt) {
      state.assign(state0);
      return false;
    }
  }

  // Advance the state from the beginning of the cycle using the midpoint derivatives.
  TIME_CheapRK2EndStep.start();
  CALI_MARK_BEGIN("TIME_CheapRK2EndStep");
  state.timeAdvanceOnly(false);
  //this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);
  this->currentTime(t + dt);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(t + dt, dt, db, state, derivs);
  this->finalizeGhostBoundaries();
  // this->enforceBoundaries(state, derivs);
  TIME_CheapRK2EndStep.stop();
  CALI_MARK_END("TIME_CheapRK2EndStep");

  // Apply any physics specific finalizations.
  TIME_CheapRK2Finalize.start();
  CALI_MARK_BEGIN("TIME_CheapRK2Finalize");
  this->postStepFinalize(t + dt, dt, state, derivs);
  TIME_CheapRK2Finalize.stop();
  CALI_MARK_END("TIME_CheapRK2Finalize");

  // Enforce boundaries.
  TIME_CheapRK2EnforceBound.start();
  CALI_MARK_BEGIN("TIME_CheapRK2EnforceBound");
  this->enforceBoundaries(state, derivs);
  TIME_CheapRK2EnforceBound.stop();
  CALI_MARK_END("TIME_CheapRK2EnforceBound");

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt);
  TIME_CheapRK2.stop();
  CALI_MARK_END("TIME_CheapRK2");

  return true;
}

}
