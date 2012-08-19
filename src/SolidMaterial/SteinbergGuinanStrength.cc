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

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::abs;
using std::min;
using std::max;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
SteinbergGuinanStrength<Dimension, Constants>::
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
  StrengthModel(),
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
template<typename Dimension, typename Constants>
SteinbergGuinanStrength<Dimension, Constants>::
~SteinbergGuinanStrength() {
}

//------------------------------------------------------------------------------
// Compute the shear modulus.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
shearModulus(const double density,
             const double specificThermalEnergy,
             const double pressure) const {

  if ((mG0 > 0.0) &&
      (mA == 0.0) &&
      (mB == 0.0) &&
      (mbeta == 0.0) &&
      (mnhard == 0.0)) {
    return mG0;

  } else {
    const double eta = mEOSPtr->boundedEta(density);
    CHECK(distinctlyGreaterThan(eta, 0.0));
    const double result = mG0*max(1.0e-10,
                                  meltAttenuation(density, specificThermalEnergy)*(1.0 + 
                                                                                   mA*max(0.0, pressure)/FastMath::CubeRootHalley2(eta) -
                                                                                   mB*(max(0.0, computeTemperature(density, specificThermalEnergy) - 300.0))));
    CHECK(distinctlyGreaterThan(result, 0.0));
    return result;
  }

}

//------------------------------------------------------------------------------
// Compute the yield strength.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
yieldStrength(const double density,
              const double specificThermalEnergy,
              const double pressure,
              const double plasticStrain,
              const double plasticStrainRate) const {

  if ((mG0 > 0.0) &&
      (mA == 0.0) &&
      (mB == 0.0) &&
      (mbeta == 0.0) &&
      (mnhard == 0.0)) {
    return mY0;

  } else {
    const double eta = mEOSPtr->boundedEta(density);
    CHECK(distinctlyGreaterThan(eta, 0.0));
    const double Yhard = min(mYmax, mY0*pow(1.0 + mbeta*(plasticStrain + mgamma0), mnhard));
    const double result = Yhard*shearModulus(density, specificThermalEnergy, pressure)/mG0;
    return result;
  }

}

//------------------------------------------------------------------------------
// Compute the full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
soundSpeed(const double density,
           const double specificThermalEnergy,
           const double pressure,
           const double fluidSoundSpeed) const {
  REQUIRE(distinctlyGreaterThan(density, 0.0));
  const double cs2 = fluidSoundSpeed*fluidSoundSpeed + 
    std::abs(4.0/3.0 * shearModulus(density, specificThermalEnergy, pressure)) / density;
  ENSURE(cs2 > 0.0);
  return std::sqrt(cs2);
}

//------------------------------------------------------------------------------
// Compute the melt attenuation.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
meltAttenuation(const double density, const double specificThermalEnergy) const {
  const double tiny = 1.0e-5;
  const double rho0 = mEOSPtr->referenceDensity();
  const double mu = mEOSPtr->boundedEta(density) - 1.0;
  CHECK(rho0 > 0.0);
  CHECK(mu >= -1.0);
  const double emelt = mMeltEnergyFit(mu)/rho0;
  CHECK(fuzzyGreaterThanOrEqual(emelt, 0.0));

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
template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
computeTemperature(const double density, const double specificThermalEnergy) const {
  const double Cv = mEOSPtr->specificHeat(density, specificThermalEnergy);
  const double rho0 = mEOSPtr->referenceDensity();
  const double mu = mEOSPtr->boundedEta(density) - 1.0;
  CHECK(Cv > 0.0);
  CHECK(rho0 > 0.0);
  CHECK(mu >= -1.0);
  return max(0.0, (max(0.0, specificThermalEnergy - mColdEnergyFit(mu))/rho0)/Cv);
}

//------------------------------------------------------------------------------
// Access the strength parameters.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
G0() const {
  return mG0;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
A() const {
  return mA;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
B() const {
  return mB;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
Y0() const {
  return mY0;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
Ymax() const {
  return mYmax;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
Yp() const {
  return mYp;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
beta() const {
  return mbeta;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
gamma0() const {
  return mgamma0;
}

template<typename Dimension, typename Constants>
double
SteinbergGuinanStrength<Dimension, Constants>::
nhard() const {
  return mnhard;
}

template<typename Dimension, typename Constants>
const NinthOrderPolynomialFit&
SteinbergGuinanStrength<Dimension, Constants>::
coldEnergyFit() const {
  return mColdEnergyFit;
}

template<typename Dimension, typename Constants>
const NinthOrderPolynomialFit&
SteinbergGuinanStrength<Dimension, Constants>::
meltEnergyFit() const {
  return mMeltEnergyFit;
}

}
}

