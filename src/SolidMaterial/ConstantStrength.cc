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
                 const double Y0,
                 const double muD,
                 const double YD):
  StrengthModel<Dimension>(),
  mmu0(mu0),
  mY0(Y0),
  mmuD(muD),
  mYD(YD),
  mEOSptr(nullptr) {
}

//------------------------------------------------------------------------------
// Constructor with SolidEquationOfState.
//------------------------------------------------------------------------------
template<typename Dimension>
ConstantStrength<Dimension>::
ConstantStrength(const double mu0,
                 const double Y0,
                 const SolidEquationOfState<Dimension>& eos,
                 const double muD,
                 const double YD):
  StrengthModel<Dimension>(),
  mmu0(mu0),
  mY0(Y0),
  mmuD(muD),
  mYD(YD),
  mEOSptr(&eos) {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ConstantStrength<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& density,
             const Field<Dimension, Scalar>& /*specificThermalEnergy*/,
             const Field<Dimension, Scalar>& /*pressure*/,
             const Field<Dimension, SymTensor>& damage) const {
  if (mEOSptr != nullptr and
      density/(mEOSptr->referenceDensity()) < mEOSptr->etamin()) {
    shearModulus = mmuD;
  } else {
    const auto n = damage.numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
      shearModulus(i) = (1.0 - Di)*mmu0 + Di*mmuD;
    }
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
  if (mEOSptr != nullptr and
      density/(mEOSptr->referenceDensity()) < mEOSptr->etamin()) {
    yieldStrength = mYD;
  } else {
    const auto n = damage.numInternalElements();
#pragma omp parallel for
    for (auto i = 0u; i < n; ++i) {
      const auto Di = std::max(0.0, std::min(1.0, damage(i).eigenValues().maxElement()));
      yieldStrength(i) = (1.0 - Di)*mY0 + Di*mYD;
    }
  }
}

//------------------------------------------------------------------------------
// Return the undamaged shear modulus
//------------------------------------------------------------------------------
template<typename Dimension>
double
ConstantStrength<Dimension>::
mu0() const {
  return mmu0;
}

//------------------------------------------------------------------------------
// Return the undamaged yield strength
//------------------------------------------------------------------------------
template<typename Dimension>
double
ConstantStrength<Dimension>::
Y0() const {
  return mY0;
}

//------------------------------------------------------------------------------
// Return the damaged shear modulus
//------------------------------------------------------------------------------
template<typename Dimension>
double
ConstantStrength<Dimension>::
muD() const {
  return mmuD;
}

//------------------------------------------------------------------------------
// Return the damaged yield strength
//------------------------------------------------------------------------------
template<typename Dimension>
double
ConstantStrength<Dimension>::
YD() const {
  return mYD;
}

}
