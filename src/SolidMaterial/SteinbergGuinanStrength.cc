//---------------------------------Spheral++----------------------------------//
// SteinbergGuinanStrength -- Implements the Steinberg-Guinan strength model.
//   ** Need reference **
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//
#include "SteinbergGuinanStrength.hh"
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
SteinbergGuinanStrength<Dimension>::
SteinbergGuinanStrength(const SolidEquationOfState<Dimension>& eos,
                        const double G0,     
                        const double A,      
                        const double B,      
                        const double Y0,     
                        const double Ymax,
                        const double Yp,
                        const double beta,
                        const double gamma0, 
                        const double nhard,
                        const NinthOrderPolynomialFit& coldEnergyFit,
                        const NinthOrderPolynomialFit& meltEnergyFit):
  StrengthModel<Dimension>(),
  mEOSPtr(&eos),
  mG0(G0),
  mA(A),
  mB(B),
  mY0(Y0),
  mYmax(Ymax),
  mYp(Yp),
  mbeta(beta),
  mgamma0(gamma0),
  mnhard(nhard),
  mColdEnergyFit(coldEnergyFit),
  mMeltEnergyFit(meltEnergyFit) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SteinbergGuinanStrength<Dimension>::
~SteinbergGuinanStrength() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
shearModulus(FieldSpace::Field<Dimension, Scalar>& shearModulus,
             const FieldSpace::Field<Dimension, Scalar>& density,
             const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
             const FieldSpace::Field<Dimension, Scalar>& pressure) const {
  if ((mG0 > 0.0) &&
      (mA == 0.0) &&
      (mB == 0.0) &&
      (mbeta == 0.0) &&
      (mnhard == 0.0)) {
    shearModulus = mG0;

  } else {
    Field<Dimension, Scalar> T("temperature", density.nodeList());
    this->computeTemperature(T, density, specificThermalEnergy);
    for (unsigned i = 0; i != density.numInternalElements(); ++i) {
      const double eta = mEOSPtr->boundedEta(density(i));
      CHECK(distinctlyGreaterThan(eta, 0.0));
      shearModulus(i) = mG0*max(1.0e-10,
                                meltAttenuation(density(i), specificThermalEnergy(i))*(1.0 + 
                                                                                       mA*max(0.0, pressure(i))/FastMath::CubeRootHalley2(eta) -
                                                                                       mB*max(0.0, T(i))));
      CHECK(distinctlyGreaterThan(shearModulus(i), 0.0));
    }
  }
}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
yieldStrength(FieldSpace::Field<Dimension, Scalar>& yieldStrength,
              const FieldSpace::Field<Dimension, Scalar>& density,
              const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
              const FieldSpace::Field<Dimension, Scalar>& pressure,
              const FieldSpace::Field<Dimension, Scalar>& plasticStrain,
              const FieldSpace::Field<Dimension, Scalar>& plasticStrainRate) const {
  if ((mG0 > 0.0) &&
      (mA == 0.0) &&
      (mB == 0.0) &&
      (mbeta == 0.0) &&
      (mnhard == 0.0)) {
    yieldStrength = mY0;

  } else {
    this->shearModulus(yieldStrength, density, specificThermalEnergy, pressure);
    for (unsigned i = 0; i != density.numInternalElements(); ++i) {
      const double eta = mEOSPtr->boundedEta(density(i));
      if (fuzzyEqual(eta, mEOSPtr->etamin())) {
        yieldStrength(i) = 0.0;
      } else {
        CHECK(distinctlyGreaterThan(eta, 0.0));
        const double Yhard = min(mYmax, mY0*pow(1.0 + mbeta*(plasticStrain(i) + mgamma0), mnhard));
        yieldStrength(i) *= Yhard/mG0;
      }
    }
  }
}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
soundSpeed(FieldSpace::Field<Dimension, Scalar>& soundSpeed,
           const FieldSpace::Field<Dimension, Scalar>& density,
           const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy,
           const FieldSpace::Field<Dimension, Scalar>& pressure,
           const FieldSpace::Field<Dimension, Scalar>& fluidSoundSpeed) const {
  Field<Dimension, Scalar> mu("shear modulus", density.nodeList());
  this->shearModulus(mu, density, specificThermalEnergy, pressure);
  for (unsigned i = 0; i != density.numInternalElements(); ++i) {
    CHECK(density(i) > 0.0);
    const double cs2 = fluidSoundSpeed(i)*fluidSoundSpeed(i) + std::abs(4.0/3.0 * mu(i)) / density(i);
    CHECK(cs2 > 0.0);
    soundSpeed(i) = std::sqrt(cs2);
  }
}

//------------------------------------------------------------------------------
// Compute the melt attenuation.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
meltAttenuation(const double density, const double specificThermalEnergy) const {
  const double tiny = 1.0e-5;
  const double rho0 = mEOSPtr->referenceDensity();
  const double mu = mEOSPtr->boundedEta(density) - 1.0;
  CHECK(rho0 > 0.0);
  CHECK(mu >= -1.0);
  const double emelt = mMeltEnergyFit(mu)/rho0;
  // CHECK(fuzzyGreaterThanOrEqual(emelt, 0.0));

  double result;
  if (fuzzyEqual(emelt, 0.0) || specificThermalEnergy < 0.0) {
    result = 1.0;
  } else {
    if (specificThermalEnergy > (1.0 - tiny)*emelt) {
      result = 0.0;
    } else {
      CHECK(!fuzzyEqual(specificThermalEnergy, emelt));
      result = exp(-mYp*specificThermalEnergy/(emelt - specificThermalEnergy));
    }
  }

  return result;
}

//------------------------------------------------------------------------------
// Compute the Steinberg-Guinan effective temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SteinbergGuinanStrength<Dimension>::
computeTemperature(FieldSpace::Field<Dimension, Scalar>& temperature,
                   const FieldSpace::Field<Dimension, Scalar>& density,
                   const FieldSpace::Field<Dimension, Scalar>& specificThermalEnergy) const {
  Field<Dimension, Scalar> cV("specific heat", density.nodeList());
  mEOSPtr->setSpecificHeat(cV, density, specificThermalEnergy);
  const double rho0 = mEOSPtr->referenceDensity();
  CHECK(rho0 > 0.0);
  for (unsigned i = 0; i != density.numInternalElements(); ++i) {
    const double mu = mEOSPtr->boundedEta(density(i)) - 1.0;
    CHECK(mu >= -1.0);
    CHECK(cV(i) > 0.0);
    const double emelt = mMeltEnergyFit(mu)/rho0;
    const double eps = max(0.0, min(emelt, specificThermalEnergy(i)));
    temperature(i) = max(0.0, eps - mColdEnergyFit(mu))*rho0/cV(i);
  }
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
G0() const {
  return mG0;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
A() const {
  return mA;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
B() const {
  return mB;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
Y0() const {
  return mY0;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
Ymax() const {
  return mYmax;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
Yp() const {
  return mYp;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
beta() const {
  return mbeta;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
gamma0() const {
  return mgamma0;
}

template<typename Dimension>
double
SteinbergGuinanStrength<Dimension>::
nhard() const {
  return mnhard;
}

template<typename Dimension>
const NinthOrderPolynomialFit&
SteinbergGuinanStrength<Dimension>::
coldEnergyFit() const {
  return mColdEnergyFit;
}

template<typename Dimension>
const NinthOrderPolynomialFit&
SteinbergGuinanStrength<Dimension>::
meltEnergyFit() const {
  return mMeltEnergyFit;
}

}
}

