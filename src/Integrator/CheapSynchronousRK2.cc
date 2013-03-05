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

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Construct the state and derivatives, and initalize the integrator.
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  state.timeAdvanceOnly(true);
  this->initialize(state, derivs);

  // Determine the minimum timestep across all packages.
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);

  // Copy the beginning of step state.
  State<Dimension> state0(state);
  state0.copyState();
  state0.timeAdvanceOnly(false);

  // Trial advance the state to the mid timestep point.
  const double hdt = 0.5*dt;
  state.update(derivs, hdt, t, hdt);
  this->currentTime(t + hdt);

  // Enforce Boundary conditions.
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
                                  
  // Do any physics specific stuff relating to the fact the state was just updated.
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Loop over the physics packages and perform any necessary initializations.
  this->preStepInitialize(t + hdt, hdt, state, derivs);

  // Zero out the stored derivatives.
  derivs.Zero();

  // Now loop over the packages and evaluate the derivatives at the
  // midpoint.
  this->evaluateDerivatives(t + hdt, hdt, db, state, derivs);
  this->finalizeDerivatives(t + hdt, hdt, db, state, derivs);

  // Advance the state from the beginning of the cycle using the midpoint 
  // derivatives.
  state.timeAdvanceOnly(false);
  this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);
  this->currentTime(t + dt);

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

