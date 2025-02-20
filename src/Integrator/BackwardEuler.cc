//---------------------------------Spheral++----------------------------------//
// BackwardEuler -- Advance the set of Physics packages in time using first
// order Backward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Mon Oct 21 14:32:05 PDT 2024
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "BackwardEuler.hh"
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
// Constructor
//------------------------------------------------------------------------------
template<typename Dimension>
BackwardEuler<Dimension>::
BackwardEuler(DataBase<Dimension>& dataBase,
              const vector<Physics<Dimension>*> physicsPackages,
              const Scalar beta,
              const Scalar tol,
              const size_t maxIterations):
  ImplicitIntegrator<Dimension>(dataBase, physicsPackages, tol),
  mBeta(beta),
  mMaxIters(maxIterations) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
BackwardEuler<Dimension>::~BackwardEuler() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
BackwardEuler<Dimension>&
BackwardEuler<Dimension>::
operator=(const BackwardEuler<Dimension>& rhs) {
  if (this != &rhs) {
    ImplicitIntegrator<Dimension>::operator=(rhs);
    mBeta = rhs.mBeta;
    mMaxIters = rhs.mMaxIters;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
BackwardEuler<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  // Get the current time and data base.
  Scalar t = this->currentTime();
  const Scalar tol = this->convergenceTolerance();
  DataBase<Dimension>& db = this->accessDataBase();

  // Initalize the integrator.
  this->preStepInitialize(state, derivs);

  // Determine the minimum timestep across all packages.
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);

  // Evaluate the beginning of step derivatives.
  this->initializeDerivatives(t, dt, state, derivs);
  derivs.Zero();
  this->evaluateDerivatives(t, dt, db, state, derivs);
  this->finalizeDerivatives(t, dt, db, state, derivs);

  // Make a copy of the initial state and derivatives
  State<Dimension> state0(state);
  StateDerivatives<Dimension> derivs0(derivs);
  state0.copyState();
  derivs0.copyState();
  CHECK(state0 == state);

  // Initial Forward Euler prediction for the end of timestep state
  state.update(derivs, dt, t, dt);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt, dt, db, state, derivs);

  // Iterate on the new state
  auto residual = 2.0*tol;
  auto k = 0u;
  while (k++ < mMaxIters and residual > tol) {
    
    // Make a copy of the current end of step state
    State<Dimension> state1(state);
    state1.copyState();
    CHECK(state1 == state);

    // Compute the derivatives for the current trial end of step state
    this->initializeDerivatives(t, dt, state, derivs);
    derivs.Zero();
    this->evaluateDerivatives(t, dt, db, state, derivs);
    this->finalizeDerivatives(t, dt, db, state, derivs);
    
    // Trial advance the state to the end of the timestep.
    state.assign(state0);
    state.update(derivs, mBeta * dt, t, dt);
    if (mBeta < 1.0) state.update(derivs0, (1.0 - mBeta) * dt, t, dt);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);

    // Check the current residual
    residual = this->computeResiduals(state, state1);
  }

  if (this->verbose() and Process::getRank() == 0) {
    cerr << "    BackwardEuler: " << (k - 1u) << "/" << mMaxIters << " for residual " << residual << endl;
  }

  if (residual > tol) {
    if (Process::getRank() == 0u) cerr << "BackwardEuler step WARNING: failed to converge on new state" << endl;
    this->currentTime(t);
    state.assign(state0);
    derivs.assign(derivs0);
    return false;
  }

  // Apply any physics specific finalizations.
  this->postStepFinalize(t + dt, dt, state, derivs);

  // Enforce boundaries.
  this->enforceBoundaries(state, derivs);

  // Set the new current time and last time step.
  this->currentTime(t + dt);
  this->currentCycle(this->currentCycle() + 1);
  this->lastDt(dt);

  return true;
}

}
