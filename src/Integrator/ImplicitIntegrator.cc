//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrator -- Abstract base class for implicit in time integrators.
//
// Created by JMO, Fri Oct 18 16:34:11 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/ImplicitIntegrator.hh"
#include "Distributed/allReduce.hh"

namespace Spheral {

using std::cout;
using std::endl;

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>::
ImplicitIntegrator(DataBase<Dimension>& dataBase,
                   const Scalar tol):
  Integrator<Dimension>(dataBase),
  mTol(tol) {
}

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
  return *this;
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
  auto result = ResidualType(0.0, "default");
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
