//---------------------------------Spheral++----------------------------------//
// MurnahanEquationOfState
//
//   P(rho) = 1/(nK) * (eta^n - 1)
//   eta = rho/rho0
//
// Created by JMO, Mon Jun  6 13:53:50 PDT 2005
//----------------------------------------------------------------------------//
#include "MurnahanEquationOfState.hh"
#include "Field/Field.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "Infrastructure/SpheralError.hh"

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
MurnahanEquationOfState<Dimension, Constants>::
MurnahanEquationOfState(const double referenceDensity,
                        const double etamin,
                        const double etamax,
                        const double n,
                        const double K,
                        const double atomicWeight,
                        const double externalPressure,
                        const double minimumPressure,
                        const double maximumPressure):
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  minimumPressure,
                                  maximumPressure),
  mn(n),
  mK(K),
  mAtomicWeight(atomicWeight),
  mExternalPressure(externalPressure),
  mCv(3.0 * referenceDensity * Constants::MolarGasConstant / atomicWeight),
  mnKi(1.0/(n*K)) {
  REQUIRE(distinctlyGreaterThan(n, 0.0));
  REQUIRE(distinctlyGreaterThan(K, 0.0));
  REQUIRE(distinctlyGreaterThan(mAtomicWeight, 0.0));
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
MurnahanEquationOfState<Dimension, Constants>::
~MurnahanEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  throw SpheralError("MurnahanEquationOfState::gamma UNIMPLEMENTED.");
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  CHECK(distinctlyGreaterThan(eta, 0.0));
  return max(this->minimumPressure(),
             min(this->maximumPressure(),
                 mnKi*(pow(eta, mn) - 1.0)));
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// This is a *hokey* definition -- have to do better if we ever really care
// about the temperature.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
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
MurnahanEquationOfState<Dimension, Constants>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  throw SpheralError("MurnahanEquationOfState::gamma UNIMPLEMENTED.");
  return 0.0;
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension, Constants>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return massDensity * computeDPDrho(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho).
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
MurnahanEquationOfState<Dimension, Constants>::
computeDPDrho(const Scalar massDensity,
              const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  if (fuzzyEqual(eta, this->etamin()) || 
      fuzzyEqual(eta, this->etamax())) return 0.0;
  const double result = pow(eta, mn - 1)/(mK*this->referenceDensity());
  ENSURE(result >= 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
bool
MurnahanEquationOfState<Dimension, Constants>::valid() const {
  return (SolidEquationOfState<Dimension>::valid() && 
          mn > 0.0 &&
          mK > 0.0 &&
          mAtomicWeight > 0.0 &&
          mCv > 0.0);
}

}
}

