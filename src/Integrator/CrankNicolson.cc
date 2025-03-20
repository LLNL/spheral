//---------------------------------Spheral++----------------------------------//
// CrankNicolson -- Advance the set of Physics packages in time using first
// order Backward Euler method.  All packages are advanced at one
// timestep simultaneously each step, i.e., synchronously.
//
// Created by JMO, Mon Oct 21 14:32:05 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/CrankNicolson.hh"
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
CrankNicolson<Dimension>::
CrankNicolson(DataBase<Dimension>& dataBase,
              const vector<Physics<Dimension>*> physicsPackages,
              const Scalar alpha,
              const Scalar tol,
              const size_t maxIterations):
  ImplicitIntegrator<Dimension>(dataBase, physicsPackages, tol),
  mAlpha(alpha),
  mMaxIterations(maxIterations),
  mNumExplicitSteps(0u),
  mNumImplicitSteps(0u) {
  this->allowDtCheck(true);
}

//------------------------------------------------------------------------------
// Take a step.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
CrankNicolson<Dimension>::
step(typename Dimension::Scalar maxTime,
     State<Dimension>& state,
     StateDerivatives<Dimension>& derivs) {

  // Get the current time and data base.
  const auto t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();
  const auto tol = this->convergenceTolerance();

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
  State<Dimension> state0(state), state1(state);
  StateDerivatives<Dimension> derivs0(derivs);
  state0.copyState();
  state1.copyState();
  derivs0.copyState();
  CHECK(state0 == state);

  // If we have not yet accrued enough previous step information to make a
  // prediction about the next state, just advance explicitly
  if (mDtMultiplier < 1.0) {

    //..........................................................................
    // Do a standard RK2 step
    ++mNumExplicitSteps;
    const auto hdt = 0.5*dt;
    
    // Trial advance the state to the mid timestep point.
    state.update(derivs, hdt, t, hdt);
    this->currentTime(t + hdt);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + hdt, hdt, db, state, derivs);

    // Evaluate the derivatives at the trial midpoint conditions.
    this->initializeDerivatives(t + hdt, hdt, state, derivs);
    derivs.Zero();
    this->evaluateDerivatives(t + hdt, hdt, db, state, derivs);
    this->finalizeDerivatives(t + hdt, hdt, db, state, derivs);

    // Advance the state from the beginning of the cycle using the midpoint 
    // derivatives.
    state.assign(state0);
    state.update(derivs, dt, t, dt);
    this->currentTime(t + dt);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);

  } else {

    // Initial Forward Euler prediction for the end of timestep state
    state.update(derivs, dt, t, dt, false);
    this->applyGhostBoundaries(state, derivs);
    this->finalizeGhostBoundaries();
    this->postStateUpdate(t + dt, dt, db, state, derivs);

    // Initial independent variable vector
    vector<double> solution1, solution;
    // state0.serializeIndependentData(solution0);
    // const auto n = solution0.size();

    // Iterate!
    auto done = false;
    size_t iterations = 0u;
    while (iterations++ < mMaxIterations and not done) {

      // Last pass on new state
      state1.assign(state);
      state1.serializeIndependentData(solution1);
      const auto n = solution1.size();

      // Derivatives at the last prediction
      this->initializeDerivatives(t + dt, dt, state, derivs);
      derivs.Zero();
      this->evaluateDerivatives(t + dt, dt, db, state, derivs);
      this->finalizeDerivatives(t + dt, dt, db, state, derivs);

      // Estimate new n+1 solution
      state.assign(state0);
      state.update(derivs0, 0.5*dt, t,          0.5*dt);
      state.update(derivs,  0.5*dt, t + 0.5*dt, 0.5*dt);

      // Are we blending old and new solutions?
      if (mAlpha > 0.0) {
        state.serializeIndependentData(solution);
        CHECK(solution.size() == n);
        for (auto i = 0u; i < n; ++i) solution[i] = mAlpha*solution1[i] + (1.0 - mAlpha)*solution[i];
        state.deserializeIndependentData(solution);
        state.update(derivs0, 0.5*dt, t,          0.5*dt, false);
        state.update(derivs,  0.5*dt, t + 0.5*dt, 0.5*dt, false);
      }

      // Finish state update
      this->applyGhostBoundaries(state, derivs);
      this->finalizeGhostBoundaries();
      this->postStateUpdate(t + dt, dt, db, state, derivs);

      // Compare for convergence
      const auto maxResidual = this->computeResiduals(state, state1);
      done = maxResidual < tol;
      // cerr << "=============> CrankNicolson: " << iterations << "/" << mMaxIterations << " : " << maxResidual << "/" << tol << endl;
    }

    // Did we succeed?
    if (iterations >= mMaxIterations) {
      state.assign(state0);
      return false;
    }
  }
  ++mNumImplicitSteps;
  
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
