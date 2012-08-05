//---------------------------------Spheral++----------------------------------//
// LinearPolynomialEquationOfState -- An equation of state approximated by a
// linear polynomial, i.e.:
//
//   P(rho, e) = A0 + A1*mu + a2*mu^2 + a3*mu^3 + (B0 + B1*mu + B2*mu^2)*e
//   mu = rho/rho0 - 1.0
//
// Created by JMO, Thu May  5 16:07:36 PDT 2005
//----------------------------------------------------------------------------//
#include <iostream>
using namespace std;

#include "LinearPolynomialEquationOfState.hh"
#include "Field/Field.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "DBC.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::min;
using std::max;
using std::abs;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Construct with the given coefficients.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
LinearPolynomialEquationOfState<Dimension, Constants>::
LinearPolynomialEquationOfState(const double referenceDensity,
                                const double etamin,
                                const double etamax,
                                const double a0,
                                const double a1,
                                const double a2,
                                const double a3,
                                const double b0,
                                const double b1,
                                const double b2,
                                const double atomicWeight,
                                const double externalPressure,
                                const double minimumPressure,
                                const double maximumPressure):
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  minimumPressure,
                                  maximumPressure),
  mA0(a0),
  mA1(a1),
  mA2(a2),
  mA3(a3),
  mB0(b0),
  mB1(b1),
  mB2(b2),
  mAtomicWeight(atomicWeight),
  mCv(3.0 * referenceDensity * Constants::MolarGasConstant / atomicWeight),
  mGamma(mB0 + 1.0),
  mExternalPressure(externalPressure) {
  REQUIRE(distinctlyGreaterThan(mAtomicWeight, 0.0));
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
LinearPolynomialEquationOfState<Dimension, Constants>::
~LinearPolynomialEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
LinearPolynomialEquationOfState<Dimension, Constants>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (int i = 0; i != Pressure.size(); ++i) {
    Pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
LinearPolynomialEquationOfState<Dimension, Constants>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (int i = 0; i != temperature.size(); ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
LinearPolynomialEquationOfState<Dimension, Constants>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(valid());
  for (int i = 0; i != specificThermalEnergy.size(); ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
LinearPolynomialEquationOfState<Dimension, Constants>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(valid());
  specificHeat = mCv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
LinearPolynomialEquationOfState<Dimension, Constants>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (int i = 0; i != soundSpeed.size(); ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
LinearPolynomialEquationOfState<Dimension, Constants>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  gamma = mGamma;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a 
// polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
LinearPolynomialEquationOfState<Dimension, Constants>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (int i = 0; i != bulkModulus.size(); ++i) {
    bulkModulus(i)=this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension, Constants>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  return max(this->minimumPressure(),
             min(this->maximumPressure(),
                 mA0 + mA1*mu + mA2*mu*mu + mA3*mu*mu*mu +
                 (mB0 + mB1*mu + mB2*mu*mu)*specificThermalEnergy - mExternalPressure));
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// This is a *hokey* definition -- have to do better if we ever really care
// about the temperature.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension, Constants>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return mCv*specificThermalEnergy + 300;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension, Constants>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  REQUIRE(valid());
  return (temperature - 300.0)/mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension, Constants>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  REQUIRE(valid());
  return mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension, Constants>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double c2 = computeDPDrho(massDensity, specificThermalEnergy);
  ENSURE(c2 >= 0.0);
  return sqrt(c2);
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension, Constants>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
LinearPolynomialEquationOfState<Dimension, Constants>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return massDensity * computeDPDrho(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho).
// 
// This turns out to be 
// \partial P       \partial P   |          P     \partial P   |
// ------------   = -------------|      + ------  -------------|
// \partial \rho    \partial \rho|_\eps   \rho^2  \partial \eps|_\rho
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
LinearPolynomialEquationOfState<Dimension, Constants>::
computeDPDrho(const Scalar massDensity,
              const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  const double mu = eta - 1.0;
  const double rho0 = this->referenceDensity();
  const double rho = rho0*eta;
  const double dPdrho_eps = std::abs(mA1 + mA2*mu + mA3*mu*mu +
                                     (mB1 + mB2*mu)*specificThermalEnergy)/rho0;
  const double Prho2 = this->pressure(massDensity, specificThermalEnergy)/(rho*rho);
  const double dPdeps_rho = mB0 + mB1*mu + mB2*mu*mu;
  const double result = dPdrho_eps + Prho2*dPdeps_rho;

  ENSURE(result >= 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
bool
LinearPolynomialEquationOfState<Dimension, Constants>::valid() const {
  return (SolidEquationOfState<Dimension>::valid() && 
          mAtomicWeight > 0.0 &&
          mCv > 0.0 &&
          mGamma > 0.0);
}

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Material/PhysicalConstants.hh"
#include "Material/MKSUnits.hh"
#include "Material/CGSUnits.hh"
#include "Material/SolarUnits.hh"
namespace Spheral {
  namespace SolidMaterial {
    template class LinearPolynomialEquationOfState<Dim<1>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<2>, Material::PhysicalConstants<Material::MKSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<3>, Material::PhysicalConstants<Material::MKSUnits> >;

    template class LinearPolynomialEquationOfState<Dim<1>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<2>, Material::PhysicalConstants<Material::CGSUnits> >;
    template class LinearPolynomialEquationOfState<Dim<3>, Material::PhysicalConstants<Material::CGSUnits> >;

    template class LinearPolynomialEquationOfState<Dim<1>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class LinearPolynomialEquationOfState<Dim<2>, Material::PhysicalConstants<Material::SolarUnits> >;
    template class LinearPolynomialEquationOfState<Dim<3>, Material::PhysicalConstants<Material::SolarUnits> >;
  }
}
