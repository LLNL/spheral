//---------------------------------Spheral++----------------------------------//
// BackwardEuler -- Advance the set of Physics packages in time using first
// order Backward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Mon Oct 21 14:32:05 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/BackwardEuler.hh"
#include "Integrator/ImplicitIntegrationVectorOperator.hh"
#include "Solvers/KINSOL.hh"

#include "DataOutput/Restart.hh"
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
              const Scalar ftol,
              const Scalar steptol,
              const size_t maxIterations):
  ImplicitIntegrator<Dimension>(dataBase, physicsPackages, ftol),
  mBeta(beta),
  mftol(ftol),
  msteptol(steptol),
  mMaxIters(maxIterations) {
  this->allowDtCheck(true);
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

  // Build a solver
  KINSOL solver;
  solver.fnormtol(mftol);
  solver.scsteptol(msteptol);
  solver.numMaxIters(mMaxIters);

  // Build the VectorOperator
  ImplicitIntegrationVectorOperator op(t, dt, mBeta, state0, derivs0, state, derivs, *this);

  // Iterate on the new state
  std::vector<double> solution;
  state.serializeIndependentData(solution);
  auto numIters = solver.solve(op, solution);
  if (numIters == size_t(mMaxIters)) {
    state.assign(state0);
    return false;
  }

  // Unpack the solution into the final state.
  // state.deserializeIndependentData(solution);
  state.update(derivs, dt, t, dt, false);
  this->applyGhostBoundaries(state, derivs);
  this->finalizeGhostBoundaries();
  this->postStateUpdate(t + dt, dt, db, state, derivs);

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
