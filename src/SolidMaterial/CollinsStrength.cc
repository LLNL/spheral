//---------------------------------Spheral++----------------------------------//
// CollinsYieldStrength -- Implements a pressure dependent yield strength
// model appropriate for geological materials.
//
// Since this is solely a yield strength model it takes another StrengthModel
// as an argument to compute the shear modulus.  Perhaps at some point we should
// just split up the ideas of what provides shear modulus and yield strength?
// 
//    See Collins, Melosh, Ivanov, 2004 Appendix, MAPS
//
// Created by JMO, Thu Jan 14 16:40:21 PST 2016
//     Based on python implementation by Megan Syal and Cody Raskin
//----------------------------------------------------------------------------//
#include "CollinsStrength.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"
#include "SolidEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {

using std::abs;
using std::min;
using std::max;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CollinsStrength<Dimension>::
CollinsStrength(const StrengthModel<Dimension>& shearModulusModel,
                const double mui,
                const double Y0,     
                const double Ym):
  StrengthModel<Dimension>(),
  mShearModulusModel(shearModulusModel),
  mmui(mui),
  mY0(Y0),
  mYm(Ym) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
CollinsStrength<Dimension>::
~CollinsStrength() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CollinsStrength<Dimension>::
shearModulus(Field<Dimension, Scalar>& shearModulus,
             const Field<Dimension, Scalar>& density,
             const Field<Dimension, Scalar>& specificThermalEnergy,
             const Field<Dimension, Scalar>& pressure) const {
  mShearModulusModel.shearModulus(shearModulus, density, specificThermalEnergy, pressure);
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CollinsStrength<Dimension>::
yieldStrength(Field<Dimension, Scalar>& yieldStrength,
              const Field<Dimension, Scalar>& density,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/,
              const Field<Dimension, Scalar>& pressure,
              const Field<Dimension, Scalar>& /*plasticStrain*/,
              const Field<Dimension, Scalar>& /*plasticStrainRate*/) const {
  const unsigned n = density.numInternalElements();
  const Scalar YdiffInv = safeInvVar(mYm - mY0);
  for (unsigned i = 0; i != n; ++i) {
    yieldStrength(i) = mY0 + mmui*pressure(i)/(1.0 + mmui*pressure(i)*YdiffInv);
  }
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
CollinsStrength<Dimension>::
soundSpeed(Field<Dimension, Scalar>& soundSpeed,
           const Field<Dimension, Scalar>& density,
           const Field<Dimension, Scalar>& specificThermalEnergy,
           const Field<Dimension, Scalar>& pressure,
           const Field<Dimension, Scalar>& fluidSoundSpeed) const {
  mShearModulusModel.soundSpeed(soundSpeed, density, specificThermalEnergy, pressure, fluidSoundSpeed);
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
const StrengthModel<Dimension>&
CollinsStrength<Dimension>::
shearModulusModel() const {
  return mShearModulusModel;
}

template<typename Dimension>
double
CollinsStrength<Dimension>::
mui() const {
  return mmui;
}

template<typename Dimension>
double
CollinsStrength<Dimension>::
Y0() const {
  return mY0;
}

template<typename Dimension>
double
CollinsStrength<Dimension>::
Ym() const {
  return mYm;
}

}
