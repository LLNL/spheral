//---------------------------------Spheral++----------------------------------//
// JohnsonCookStrength -- Implements the Johnson-Cook strength model.
//   ** Need reference **
//
// Created by JMO, Tue Jul 28 13:39:50 PDT 2015
//----------------------------------------------------------------------------//
#include "JohnsonCookStrength.hh"
#include "Utilities/FastMath.hh"
#include "Utilities/DBC.hh"
#include "SolidEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::abs;
using std::min;
using std::max;

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookStrength<Dimension>::
JohnsonCookStrength(const SolidEquationOfState<Dimension>& eos,
                    const StrengthModel<Dimension>& shearModulusModel,
                    const double A,
                    const double B,
                    const double C,
                    const double C4,
                    const double m,
                    const double nhard,
                    const double epsdot0,
                    const double epsdotmin,
                    const double Tmelt,
                    const double Troom):
  mEOSPtr(&eos),
  mShearModulusModelPtr(&shearModulusModel),
  mA(A),
  mB(B),
  mC(C),
  mC4(C4),
  mm(m),
  mnhard(nhard),
  mEpsdot0(epsdot0),
  mEpsdotmin(epsdotmin),
  mTmelt(Tmelt),
  mTroom(Troom) {
  VERIFY2(mEpsdot0 > 0.0, "JohnsonCookStrength ERROR: reference strain-rate must be greater than zero.");
  VERIFY2(mTmelt > mTroom, "JohnsonCookStrength ERROR: Tmelt must be greater than or equal Troom.");
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
JohnsonCookStrength<Dimension>::
~JohnsonCookStrength() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookStrength<Dimension>::
shearModulus(FieldSpace::Field<Dimension, Scalar>& shearModulus,
             const FieldSpace::Field<Dimension, Scalar>& density,
             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
             const FieldSpace::Field<Dimension, Scalar>& pressure) const {
  mShearModulusModelPtr->shearModulus(shearModulus, density, specificThermalEnergy, pressure);
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookStrength<Dimension>::
yieldStrength(FieldSpace::Field<Dimension, Scalar>& yieldStrength,
              const FieldSpace::Field<Dimension, Scalar>& density,
              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
              const FieldSpace::Field<Dimension, Scalar>& pressure,
              const FieldSpace::Field<Dimension, Scalar>& plasticStrain,
              const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate) const {
  Field<Dimension, Scalar> T("temperature", yieldStrength.nodeList());
  mEOSPtr->setTemperature(T, density, specificThermalEnergy);
  for (unsigned i = 0; i != yieldStrength.numInternalElements(); ++i) {
    const double Tstar = max(0.0, T(i) - mTroom)/(mTmelt - mTroom);
    yieldStrength(i) = 
      (mA + mB*pow(plasticStrain(i), mnhard))*
      (1.0 + mC*log(max(mEpsdotmin, plasticStrainRate(i))/mEpsdot0))*
      (1.0 - pow(Tstar, mm)) +
      mC4*pressure(i);
  }
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
JohnsonCookStrength<Dimension>::
soundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
           const FieldSpace::Field<Dimension, Scalar>& density,
           const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
           const FieldSpace::Field<Dimension, Scalar>& pressure,
           const FieldSpace::Field<Dimension, Scalar>& fluidSoundSpeed) const {
  mShearModulusModelPtr->soundSpeed(soundSpeed, density, specificThermalEnergy, pressure, fluidSoundSpeed);
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
C() const {
  return mC;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
C4() const {
  return mC4;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
m() const {
  return mm;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
nhard() const {
  return mnhard;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
epsdot0() const {
  return mEpsdot0;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
epsdotmin() const {
  return mEpsdotmin;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
Tmelt() const {
  return mTmelt;
}

template<typename Dimension>
double
JohnsonCookStrength<Dimension>::
Troom() const {
  return mTroom;
}

}
}

