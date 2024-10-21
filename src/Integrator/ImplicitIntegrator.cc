//---------------------------------Spheral++----------------------------------//
// ImplicitIntegrator -- Abstract base class for implicit in time integrators.
//
// Created by JMO, Fri Oct 18 16:34:11 PDT 2024
//----------------------------------------------------------------------------//
#include "Integrator/ImplicitIntegrator.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>::
ImplicitIntegrator():
  Integrator<Dimension>() {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>::
ImplicitIntegrator(DataBase<Dimension>& dataBase):
  Integrator<Dimension>(dataBase) {
}

//------------------------------------------------------------------------------
// Construct with the given DataBase and Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
ImplicitIntegrator<Dimension>::
ImplicitIntegrator(DataBase<Dimension>& dataBase,
                   const std::vector<Physics<Dimension>*>& physicsPackages):
  Integrator<Dimension>(dataBase, physicsPackages) {
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
