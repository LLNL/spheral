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
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Initalize the integrator.
  this->preStepInitialize(state, derivs);
  // this->initializeDerivatives(t, 0.0, state, derivs);

  // Determine the minimum timestep across all packages.
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);
  const double hdt = 0.5*dt;

  // Copy the beginning of step state.
  State<Dimension> state0(state);
  state0.copyState();

  // Trial advance the state to the mid timestep point.
  state.timeAdvanceOnly(true);
  state.update(derivs, hdt, t, hdt);
  this->currentTime(t + hdt);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(t + hdt, hdt, db, state, derivs);
  this->finalizeGhostBoundaries();

  // Evaluate the derivatives at the midpoint.
  this->initializeDerivatives(t + hdt, hdt, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t + hdt, hdt, db, state, derivs);
  this->finalizeDerivatives(t + hdt, hdt, db, state, derivs);

  // Advance the state from the beginning of the cycle using the midpoint derivatives.
  state.timeAdvanceOnly(false);
  //this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);
  this->currentTime(t + dt);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(t + dt, dt, db, state, derivs);
  this->finalizeGhostBoundaries();
  // this->enforceBoundaries(state, derivs);

  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
}

}
