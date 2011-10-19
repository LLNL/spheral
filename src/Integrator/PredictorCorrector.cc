//---------------------------------Spheral++----------------------------------//
// Second order predictor corrector time integrator.
// Based on the predictor corrector scheme described in Monaghans tensile
// instability paper (Monaghan 2000, JCP 159, 290.)
//
// Created by JMO, Wed Dec  4 21:39:36 PST 2002
//----------------------------------------------------------------------------//
#include "DataOutput/Restart.hh"
#include "PredictorCorrector.hh"
#include "DataBase/DataBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Field/NodeIterators.hh"
#include "Physics/Physics.hh"

#include "DBC.hh"
#include "cdebug.hh"

#include "TAU.h"

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
PredictorCorrector<Dimension>::PredictorCorrector():
  Integrator<Dimension>() {
  cdebug << "PredictorCorrector::PredictorCorrector()" << endl;
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>::
PredictorCorrector(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
  cdebug << "PredictorCorrector::PredictorCorrector(DataBase)" << endl;
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>::
PredictorCorrector(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
  cdebug << "PredictorCorrector::PredictorCorrector(DataBase, PhysicsPackages)" << endl;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>::~PredictorCorrector() {
  cdebug << "PredictorCorrector::~PredictorCorrector()" << endl;
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
PredictorCorrector<Dimension>&
PredictorCorrector<Dimension>::
operator=(const PredictorCorrector<Dimension>& rhs) {
  cdebug << "PredictorCorrector::operator=()" << endl;
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
PredictorCorrector<Dimension>::
step(typename Dimension::Scalar maxTime) {

  // TAU timers.
  TAU_PROFILE("PredictorCorrector", "::step", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2ConstructState, "PredictorCorrector", "::step : Construct state", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2Dt, "PredictorCorrector", "::step : Set dt", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2Zero1, "PredictorCorrector", "::step : Zero derivs  1", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2EvalDerivs1, "PredictorCorrector", "::step : Eval derivs mid step", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2CopyFields, "PredictorCorrector", "::step : Copy initial state", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2PredStep, "PredictorCorrector", "::step : Predictor step", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2Boundaries1, "PredictorCorrector", "::step : Set boundaries 1", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2PostState1, "PredictorCorrector", "::step : post state update 1", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2StepInitialize1, "PredictorCorrector", "::step : Pre-step initialize 1", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2Zero2, "PredictorCorrector", "::step : Zero derivs  2", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2EvalDerivs2, "PredictorCorrector", "::step : Eval derivs end step", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2CorStep, "PredictorCorrector", "::step : Corrector step", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2Boundaries2, "PredictorCorrector", "::step : Enforce boundaries 2", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2PostState2, "PredictorCorrector", "::step : post state update 2", TAU_USER);
  TAU_PROFILE_TIMER(TimePC2Finalize, "PredictorCorrector", "::step : finalize physics", TAU_USER);

  cdebug << "PredictorCorrector::step" << endl;

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Construct the state and derivatives, and initalize the integrator.
  TAU_PROFILE_START(TimePC2ConstructState);
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  this->initialize(state, derivs);
  TAU_PROFILE_STOP(TimePC2ConstructState);

  // Determine the minimum timestep across all packages.
  TAU_PROFILE_START(TimePC2Dt);
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);
  cdebug << "PredictorCorrector::step: chose dt = " << dt << endl;
  TAU_PROFILE_STOP(TimePC2Dt);

  // Zero out the derivatives.
  TAU_PROFILE_START(TimePC2Zero1);
  derivs.Zero();
  TAU_PROFILE_STOP(TimePC2Zero1);

  // Evaluate the beginning of step derivatives.
  TAU_PROFILE_START(TimePC2EvalDerivs1);
  this->evaluateDerivatives(t, dt, db, state, derivs);
  this->finalizeDerivatives(t, dt, db, state, derivs);
  TAU_PROFILE_STOP(TimePC2EvalDerivs1);

  // Copy the beginning of step state and derivatives.
  TAU_PROFILE_START(TimePC2CopyFields);
  State<Dimension> state0(state);
  StateDerivatives<Dimension> derivs0(derivs);
  state0.copyState();
  derivs0.copyState();
  TAU_PROFILE_STOP(TimePC2CopyFields);

  // Predictor step -- trial advance to the end of step.
  TAU_PROFILE_START(TimePC2PredStep);
  state.update(derivs, dt, t, dt);
  TAU_PROFILE_STOP(TimePC2PredStep);

  // Enforce Boundary conditions.
  TAU_PROFILE_START(TimePC2Boundaries1);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  TAU_PROFILE_STOP(TimePC2Boundaries1);
                               
  // Do any physics specific stuff relating to the fact the state was just updated.
  TAU_PROFILE_START(TimePC2PostState1);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();
  TAU_PROFILE_STOP(TimePC2PostState1);

  // Loop over the physics packages and perform any necessary initializations.
  TAU_PROFILE_START(TimePC2StepInitialize1);
  this->preStepInitialize(t + dt, dt, state, derivs);
  TAU_PROFILE_STOP(TimePC2StepInitialize1);

  // Zero out the stored derivatives.
  TAU_PROFILE_START(TimePC2Zero2);
  derivs.Zero();
  TAU_PROFILE_STOP(TimePC2Zero2);

  // Now loop over the packages and evaluate the derivatives at the
  // predicted end of step state.
  TAU_PROFILE_START(TimePC2EvalDerivs2);
  this->evaluateDerivatives(t + dt, dt, db, state, derivs);
  this->finalizeDerivatives(t + dt, dt, db, state, derivs);
  TAU_PROFILE_STOP(TimePC2EvalDerivs2);

  // Apply the corrector stage.
  TAU_PROFILE_START(TimePC2CorStep);
  const Scalar hdt = 0.5*dt;
  this->copyGhostState(state, state0);
  state.assign(state0);
  state.update(derivs0, hdt, t, hdt);
  this->applyGhostBoundaries(state, derivs);
  state.update(derivs, hdt, t + hdt, hdt);
  TAU_PROFILE_STOP(TimePC2CorStep);

  // Enforce boundaries.
  TAU_PROFILE_START(TimePC2Boundaries2);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  TAU_PROFILE_STOP(TimePC2Boundaries2);

  // Do any physics specific stuff relating to the fact the state was just updated.
  TAU_PROFILE_START(TimePC2PostState2);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();
  TAU_PROFILE_STOP(TimePC2PostState2);

  // Apply any physics specific finalizations.
  TAU_PROFILE_START(TimePC2Finalize);
  this->finalize(t + dt, dt, state, derivs);
  TAU_PROFILE_STOP(TimePC2Finalize);

  // Set the new current time and last time step.
  this->currentCycle(this->currentCycle() + 1);
  this->currentTime(t + dt);
  this->lastDt(dt);
}
}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
namespace IntegratorSpace {
template class PredictorCorrector< Dim<1> >;
template class PredictorCorrector< Dim<2> >;
template class PredictorCorrector< Dim<3> >;
}
}
