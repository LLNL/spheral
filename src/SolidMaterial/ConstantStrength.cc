//---------------------------------Spheral++----------------------------------//
// ConstantStrength -- An implentation of StrengthModel returning constant
// values for the shear modulus and yield strength.
//
// Created by JMO, Mon Jul 24 10:34:40 PDT 2006
//----------------------------------------------------------------------------//
#include "ConstantStrength.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantStrength<Dimension>::
ConstantStrength(const double mu0,
                 const double Y0):
  StrengthModel<Dimension>(),
  mShearModulus0(mu0),
  mYieldStrength0(Y0),
  mEOSptr(0) {
}

//------------------------------------------------------------------------------
// Constructor with SolidEquationOfState.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantStrength<Dimension>::
ConstantStrength(const double mu0,
                 const double Y0,
                 const SolidEquationOfState<Dimension>& eos):
  StrengthModel<Dimension>(),
  mShearModulus0(mu0),
  mYieldStrength0(Y0),
  mEOSptr(&eos) {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantStrength<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& /*density*/,
             const Field<Dimension, Scalar>& /*specificThermalEnergy*/,
             const Field<Dimension, Scalar>& /*pressure*/,
             const Field<Dimension, SymTensor>& damage) const {
  const auto n = shearModulus.numInternalElements();
#pragma omp for
  for (auto i = 0u; i < n; ++i) {
    shearModulus(i) = mShearModulus0 * std::max(0.0, std::min(1.0, 1.0 - damage(i).eigenValues().maxElement()));
  }
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantStrength<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/,
              const Field<Dimension, Scalar>& /*pressure*/,
              const Field<Dimension, Scalar>& /*plasticStrain*/,
              const Field<Dimension, Scalar>& /*plasticStrainRate*/,
              const Field<Dimension, SymTensor>& damage) const {
  if (mEOSptr != 0 and
      density/(mEOSptr->referenceDensity()) < mEOSptr->etamin()) {
    yieldStrength = 0.0;
  } else {
    const auto n = yieldStrength.numInternalElements();
#pragma omp for
    for (auto i = 0u; i < n; ++i) {
      yieldStrength(i) = mYieldStrength0 * std::max(0.0, std::min(1.0, 1.0 - damage(i).eigenValues().maxElement()));
    }
  }
}

//------------------------------------------------------------------------------
// Return the input shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
double
ConstantStrength<Dimension>::
mu0() const {
  return mShearModulus0;
}

//------------------------------------------------------------------------------
// Return the input yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
double
ConstantStrength<Dimension>::
Y0() const {
  return mYieldStrength0;
}

}
