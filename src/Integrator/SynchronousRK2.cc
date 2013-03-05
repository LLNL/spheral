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
  cdebug << "SynchronousRK2::step: chose dt = " << dt << endl;

  // Zero out the derivatives.
  derivs.Zero();

  // Evaluate the beginning of step derivatives.
  const double hdt = 0.5*dt;
  this->evaluateDerivatives(t, hdt, db, state, derivs);
  this->finalizeDerivatives(t, hdt, db, state, derivs);

  // Copy the beginning of step state.
  State<Dimension> state0(state);
  state0.copyState();

  // Trial advance the state to the mid timestep point.
  state.update(derivs, hdt, t, hdt);

  // Enforce Boundary conditions (as a side effect this updates the
  // neighboring).
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
  this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs, dt, t, dt);

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

