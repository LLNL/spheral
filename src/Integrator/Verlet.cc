//---------------------------------Spheral++----------------------------------//
// Verlet -- Advance the set of Physics packages in time using the second
// order Verlet algorithm.  This method is symplectic in the absence of 
// dissipation.
//
// Based on the description in 
// Monaghan JJ. Smoothed particle hydrodynamics. Reports on progress in physics. 2005
//
// Created by JMO, Sat Aug 23 10:33:33 PDT 2014
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "Verlet.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementFieldList.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"
#include "Hydro/HydroFieldNames.hh"

#include "Utilities/DBC.hh"

namespace Spheral {
namespace IntegratorSpace {

using namespace std;

using DataBaseSpace::DataBase;
using FieldSpace::Field;
using FieldSpace::FieldList;
using PhysicsSpace::Physics;

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
Verlet<Dimension>::Verlet():
  Integrator<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
Verlet<Dimension>::
Verlet(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
Verlet<Dimension>::
Verlet(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
Verlet<Dimension>::~Verlet() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
Verlet<Dimension>&
Verlet<Dimension>::
operator=(const Verlet<Dimension>& rhs) {
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
Verlet<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  typedef typename State<Dimension>::PolicyPointer PolicyPointer;

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Initalize the integrator.
  this->preStepInitialize(state, derivs);

  // Extract velocity into its own State as a copy.
  State<Dimension> velocityState;
  FieldList<Dimension, Vector> vel = state.fields(HydroFieldNames::velocity, Vector::zero);
  const typename State<Dimension>::KeyType velKey = State<Dimension>::key(vel);
  PolicyPointer velPolicy = state.policy(velKey);
  velocityState.enroll(vel, velPolicy);
  velocityState.copyState();

  // Determine the minimum timestep across all packages.
  const Scalar dtMin = min(this->dtMin(), maxTime - t);
  const Scalar dtMax = min(this->dtMax(), maxTime - t);
  const Scalar dt0 = this->selectDt(dtMin, dtMax, state, derivs);
  const Scalar hdt0 = 0.5*dt0;

  // Evaluate the beginning of step derivatives.
  derivs.Zero();
  this->initializeDerivatives(t, hdt0, state, derivs);
  this->evaluateDerivatives(t, hdt0, db, state, derivs);
  this->finalizeDerivatives(t, hdt0, db, state, derivs);

  // Predict state to the mid-point.
  state.update(derivs, hdt0, t, dt0);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Evaluate the derivatives at the midpoint.
  this->initializeDerivatives(t + hdt0, hdt0, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t + hdt0, hdt0, db, state, derivs);
  this->finalizeDerivatives(t + hdt0, hdt0, db, state, derivs);

  // Update the independent velocity from the beginning of step to the midpoint using
  // the midpoint derivatives.
  velocityState.update(derivs, hdt0, t, dt0);

  // Copy the new mid-point velocity estimate back into the main state.  We also assign DxDt
  // in the derivatives to be the new velocity.  Note this is inconsistent with XSPH type
  // approximations!
  FieldList<Dimension, Vector> velCopy = velocityState.fields(HydroFieldNames::velocity, Vector::zero);
  FieldList<Dimension, Vector> DxDt = derivs.fields(IncrementFieldList<Dimension, Field<Dimension, Vector> >::prefix() + HydroFieldNames::position, Vector::zero);
  vel.assignFields(velCopy);
  DxDt.assignFields(velCopy);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();

  // Figure out the new time-step for the second half of our advance.
  const Scalar dt12 = this->selectDt(dtMin, dtMax, state, derivs);
  const Scalar dt1 = min(dtMax, max(dtMin, 1.0/(2.0/dt12 - 1.0/dt0)));
  const Scalar hdt1 = min(0.5*dt1, maxTime - t - hdt0);
  const Scalar dt = hdt0 + hdt1;

  // Advance the velocity to the end of step.
  velocityState.update(derivs, hdt1, t + hdt0, dt1);

  // Copy the new final velocity into the main state.
  vel.assignFields(velCopy);
  DxDt.assignFields(velCopy);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();

  // Now advance everything to the end of step.
  state.update(derivs, hdt1, t, dt1);

  // Reset the velocity back to t1 values.
  vel.assignFields(velCopy);
  DxDt.assignFields(velCopy);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);

  // Do any physics specific stuff relating to the fact the state was just updated.
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
}
}
}

