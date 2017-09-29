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
  auto  t = this->currentTime();
  auto& db = this->accessDataBase();

  // Initalize the integrator.
  this->preStepInitialize(state, derivs);

  // Copy the beginning of step positions.
  auto pos0 = state.fields(HydroFieldNames::position, Vector::zero);
  pos0.copyFields();

  // Evaluate the beginning of step derivatives.
  // The zeros here are the timestep estimate, which we're putting off 'til we do the first derivative evaluation.
  this->initializeDerivatives(t, 0.0, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t, 0.0, db, state, derivs);
  this->finalizeDerivatives(t, 0.0, db, state, derivs);

  // Determine the minimum timestep across all packages.
  const Scalar dtMin = min(this->dtMin(), maxTime - t);
  const Scalar dtMax = min(this->dtMax(), maxTime - t);
  const Scalar dt0 = this->selectDt(dtMin, dtMax, state, derivs);
  const Scalar hdt0 = 0.5*dt0;

  // Predict state at the mid-point.
  state.timeAdvanceOnly(true);
  state.update(derivs, hdt0, t, dt0);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Copy the mid-point state.
  State<Dimension> state12(state);

  // Advance the position to the end of step using the half-step velocity.
  auto vel12 = state.fields(HydroFieldNames::velocity, Vector::zero);
  pos0 += dt0*vel12;

  // Predict state at the end-point, but override the positions with our time-centered prediction.
  state.update(derivs, hdt0, t, dt0);
  {
    auto pos = state.fields(HydroFieldNames::position, Vector::zero);
    pos.assignFields(pos0);
  }
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Evaluate the derivatives at the predicted end-point.
  this->currentTime(t + dt0);
  this->initializeDerivatives(t + dt0, dt0, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t + dt0, dt0, db, state, derivs);
  this->finalizeDerivatives(t + dt0, dt0, db, state, derivs);

  // Correct the final state by the end-point derivatives.
  state.assign(state12);
  state.timeAdvanceOnly(false);
  state.update(derivs, hdt0, t + hdt0, dt0);
  {
    auto pos = state.fields(HydroFieldNames::position, Vector::zero);
    pos.assignFields(pos0);
  }
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();

  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt0, dt0, state, derivs);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt0);
  this->lastDt(dt0);
}
}
}

