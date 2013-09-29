//---------------------------------Spheral++----------------------------------//
// SynchronousRK1 -- Advance the set of Physics packages in time using first
// order Runge-Kutta -- i.e., forward Euler.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Tue Aug  1 14:44:27 PDT 2006
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "SynchronousRK1.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"

#include "Utilities/DBC.hh"

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
SynchronousRK1<Dimension>::SynchronousRK1():
  Integrator<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>::
SynchronousRK1(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>::
SynchronousRK1(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>::~SynchronousRK1() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>&
SynchronousRK1<Dimension>::
operator=(const SynchronousRK1<Dimension>& rhs) {
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
SynchronousRK1<Dimension>::
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

  // Advance the state to the end of the timestep.
  state.update(derivs, dt, t, dt);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);

  // Do any physics specific stuff relating to the fact the state was just updated.
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Apply any physics specific finalizations.
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = this->physicsPackagesBegin();
       physicsItr != this->physicsPackagesEnd();
       ++physicsItr) (*physicsItr)->finalize(t + dt, dt, db, state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
}
}
}

