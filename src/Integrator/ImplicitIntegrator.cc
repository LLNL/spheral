//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrator -- Abstract base class for implicit in time integrators.
//
// Created by JMO, Fri Oct 18 16:34:11 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/ImplicitIntegrator.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Distributed/allReduce.hh"

namespace Spheral {

using std::cout;
using std::endl;
using std::min;

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>::
ImplicitIntegrator(DataBase<Dimension>& dataBase,
                   const std::vector<Physics<Dimension>*>& physicsPackages,
                   const Scalar tol):
  Integrator<Dimension>(dataBase, physicsPackages),
  mTol(tol) {
}

//------------------------------------------------------------------------------
// Destructor
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>::
~ImplicitIntegrator() {
}

//------------------------------------------------------------------------------
// Assignment
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>&
ImplicitIntegrator<Dimension>::
operator=(const ImplicitIntegrator<Dimension>& rhs) {
  Integrator<Dimension>::operator=(rhs);
  mTol = rhs.mTol;
  return *this;
}

//------------------------------------------------------------------------------
// step
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ImplicitIntegrator<Dimension>::
step(const typename Dimension::Scalar maxTime) {

  DataBase<Dimension>& db = this->accessDataBase();
  State<Dimension> state(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  StateDerivatives<Dimension> derivs(db, this->physicsPackagesBegin(), this->physicsPackagesEnd());
  auto success = false;
  auto count = 0u;
  auto maxIterations = 10u;
  while (not success and count++ < maxIterations) {
    cerr << "=============> Current dtMultiplier = " << mDtMultiplier << endl;
    
    // Try to advance using the current timestep multiplier
    success = this->step(maxTime, state, derivs);

    // Adjust the current timestep multiplier based on whether we succeeded or not
    mDtMultiplier *= (success ?
                      1.2 :
                      0.8);
    // mDtMultiplier = min(1.0, mDtMultiplier);

    if (not success and
        this->verbose() and
        Process::getRank() == 0) {
      cerr << "ImplicitIntegrator::step did not converge with tolerance on iteration " << count << "/" << maxIterations << endl
           << "                         reducing timestep multiplier to " << mDtMultiplier << endl;
    }
  }

  // VERIFY(success);
  return success;
}

//------------------------------------------------------------------------------
// Find the maximum residual difference in the states
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
ImplicitIntegrator<Dimension>::
computeResiduals(const State<Dimension>& state1,
                 const State<Dimension>& state0) const {

  // Get the local (to this MPI rank) answer
  const auto& db = this->dataBase();
  const auto& packages = this->physicsPackages();
  const auto  tol = this->convergenceTolerance();
  auto result = ResidualType(-1.0, "DEFAULT : You should not get this");
  for (const auto* pkg: packages) {
    const auto pres = pkg->maxResidual(db, state1, state0, tol);
    if (pres.first > result.first) result = pres;
  }

  // Reduce for the global result, and optionally print out some verbose info
  const auto globalMax = allReduce(result.first, SPHERAL_OP_MAX);
  if (result.first == globalMax and this->verbose()) {
    cout << "Global residual of "
         << result.first << endl
         << result.second << endl;
  }
  cout.flush();
  return globalMax;
}

//------------------------------------------------------------------------------
// How should we query a physics package for the time step?
//------------------------------------------------------------------------------
template<typename Dimension>
typename ImplicitIntegrator<Dimension>::TimeStepType
ImplicitIntegrator<Dimension>::
dt(const Physics<Dimension>* pkg,
   const DataBase<Dimension>& dataBase,
   const State<Dimension>& state,
   const StateDerivatives<Dimension>& derivs,
   const Scalar currentTime) const {
  return pkg->dtImplicit(dataBase, state, derivs, currentTime);
}

}
