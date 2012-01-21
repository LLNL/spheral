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
SynchronousRK1<Dimension>::SynchronousRK1():
  Integrator<Dimension>() {
  cdebug << "SynchronousRK1::SynchronousRK1()" << endl;
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>::
SynchronousRK1(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
  cdebug << "SynchronousRK1::SynchronousRK1(DataBase)" << endl;
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>::
SynchronousRK1(DataBase<Dimension>& dataBase,
               const vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
  cdebug << "SynchronousRK1::SynchronousRK1(DataBase, PhysicsPackages)" << endl;
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>::~SynchronousRK1() {
  cdebug << "SynchronousRK1::~SynchronousRK1()" << endl;
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
SynchronousRK1<Dimension>&
SynchronousRK1<Dimension>::
operator=(const SynchronousRK1<Dimension>& rhs) {
  cdebug << "SynchronousRK1::operator=()" << endl;
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

  // TAU timers.
  TAU_PROFILE("SynchronousRK1", "::step", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1ConstructState, "SynchronousRK1", "::step : Construct state ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1Dt, "SynchronousRK1", "::step : Set dt", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1Zero1, "SynchronousRK1", "::step : Zero derivs  1          ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1EvalDerivs1, "SynchronousRK1", "::step : Eval derivs @ t0  ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1MidStep, "SynchronousRK1", "::step : Advance to mid-step", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1Boundaries2, "SynchronousRK1", "::step : Set boundaries 2        ", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1PostState1, "SynchronousRK1", "::step: Post state update 1", TAU_USER);
  TAU_PROFILE_TIMER(TimeRK1Finalize, "SynchronousRK1", "::step : finalize physics", TAU_USER);

  // Get the current time and data base.
  Scalar t = this->currentTime();
  DataBase<Dimension>& db = this->accessDataBase();

  // Construct the state and derivatives, and initalize the integrator.
  TAU_PROFILE_START(TimeRK1ConstructState);
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  this->initialize(state, derivs);
  TAU_PROFILE_STOP(TimeRK1ConstructState);

  // Determine the minimum timestep across all packages.
  TAU_PROFILE_START(TimeRK1Dt);
  const Scalar dt = this->selectDt(min(this->dtMin(), maxTime - t),
                                   min(this->dtMax(), maxTime - t),
                                   state,
                                   derivs);
  cdebug << "SynchronousRK1::step: chose dt = " << dt << endl;
  TAU_PROFILE_STOP(TimeRK1Dt);

  // Zero out the derivatives.
  TAU_PROFILE_START(TimeRK1Zero1);
  derivs.Zero();
  TAU_PROFILE_STOP(TimeRK1Zero1);

  // Evaluate the beginning of step derivatives.
  TAU_PROFILE_START(TimeRK1EvalDerivs1);
  this->evaluateDerivatives(t, dt, db, state, derivs);
  this->finalizeDerivatives(t, dt, db, state, derivs);
  TAU_PROFILE_STOP(TimeRK1EvalDerivs1);

  // Advance the state to the end of the timestep.
  TAU_PROFILE_START(TimeRK1MidStep);
  state.update(derivs, dt, t, dt);
  TAU_PROFILE_STOP(TimeRK1MidStep);

  // Enforce boundaries.
  TAU_PROFILE_START(TimeRK1Boundaries2);
  this->enforceBoundaries(state, derivs);
  this->applyGhostBoundaries(state, derivs);
  TAU_PROFILE_STOP(TimeRK1Boundaries2);

  // Do any physics specific stuff relating to the fact the state was just updated.
  TAU_PROFILE_START(TimeRK1PostState1);
  this->postStateUpdate(db, state, derivs);
  this->finalizeGhostBoundaries();
  TAU_PROFILE_STOP(TimeRK1PostState1);

  // Apply any physics specific finalizations.
  TAU_PROFILE_START(TimeRK1Finalize);
  for (typename Integrator<Dimension>::ConstPackageIterator physicsItr = this->physicsPackagesBegin();
       physicsItr != this->physicsPackagesEnd();
       ++physicsItr) (*physicsItr)->finalize(t + dt, dt, db, state, derivs);
  TAU_PROFILE_STOP(TimeRK1Finalize);

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
template class SynchronousRK1< Dim<1> >;
template class SynchronousRK1< Dim<2> >;
template class SynchronousRK1< Dim<3> >;
}
}
