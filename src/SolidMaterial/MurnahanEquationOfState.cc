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
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Construct with the given coefficients.
//------------------------------------------------------------------------------
template<typename Dimension>
MurnahanEquationOfState<Dimension>::
MurnahanEquationOfState(const double referenceDensity,
                        const double etamin,
                        const double etamax,
                        const double n,
                        const double K,
                        const double atomicWeight,
                        const PhysicalConstants& constants,
                        const double externalPressure,
                        const double minimumPressure,
                        const double maximumPressure,
                        const MaterialPressureMinType minPressureType):
  SolidEquationOfState<Dimension>(referenceDensity,
                                  etamin,
                                  etamax,
                                  constants,
                                  minimumPressure,
                                  maximumPressure,
                                  minPressureType),
  mn(n),
  mK(K),
  mAtomicWeight(atomicWeight),
  mExternalPressure(externalPressure),
  mCv(3.0 * constants.molarGasConstant() / atomicWeight),
  mnKi(1.0/(n*K)) {
  REQUIRE(distinctlyGreaterThan(n, 0.0));
  REQUIRE(distinctlyGreaterThan(K, 0.0));
  REQUIRE(distinctlyGreaterThan(mAtomicWeight, 0.0));
  ENSURE(valid());
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
MurnahanEquationOfState<Dimension>::
~MurnahanEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (auto i = 0u; i != Pressure.size(); ++i) {
    Pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (auto i = 0u; i != temperature.size(); ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(valid());
  for (auto i = 0u; i != specificThermalEnergy.size(); ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  REQUIRE(valid());
  specificHeat = mCv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (auto i = 0u; i != soundSpeed.size(); ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (auto i = 0u; i != gamma.size(); ++i) {
    gamma(i) = this->gamma(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (auto i = 0u; i != bulkModulus.size(); ++i) {
    bulkModulus(i)=this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
MurnahanEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    entropy(i) = pressure(massDensity(i), specificThermalEnergy(i))*safeInvVar(pow(massDensity(i), gamma(massDensity(i), specificThermalEnergy(i))));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
pressure(const Scalar massDensity,
         const Scalar /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  if (fuzzyEqual(eta, this->etamin())) return 0.0;
  CHECK(distinctlyGreaterThan(eta, 0.0));
  return this->applyPressureLimits(mnKi*(pow(eta, mn) - 1.0));
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
// This is a *hokey* definition -- have to do better if we ever really care
// about the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
temperature(const Scalar /*massDensity*/,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return specificThermalEnergy/mCv + 300;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
specificThermalEnergy(const Scalar /*massDensity*/,
                      const Scalar temperature) const {
  REQUIRE(valid());
  return (temperature - 300.0)*mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
specificHeat(const Scalar /*massDensity*/,
             const Scalar /*temperature*/) const {
  REQUIRE(valid());
  return mCv;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
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
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
gamma(const Scalar massDensity,
      const Scalar /*specificThermalEnergy*/) const {
  const double eta = this->boundedEta(massDensity),
               rho0 = this->referenceDensity(),
               rho = rho0*eta,
               nDen = rho/mAtomicWeight;
  CHECK(mCv > 0.0);
  return 1.0 + mConstants.molarGasConstant()*nDen/mCv;
}

//------------------------------------------------------------------------------
// Calculate the individual bulk modulus.  
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return massDensity * computeDPDrho(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
MurnahanEquationOfState<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return this->pressure(massDensity, specificThermalEnergy)*safeInvVar(pow(massDensity, gamma(massDensity, specificThermalEnergy)));
}

//------------------------------------------------------------------------------
// Compute (\partial P)/(\partial rho).
//------------------------------------------------------------------------------
template<typename Dimension>
double
MurnahanEquationOfState<Dimension>::
computeDPDrho(const Scalar massDensity,
              const Scalar /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  const double eta = this->boundedEta(massDensity);
  if (fuzzyEqual(eta, this->etamin()) || 
      fuzzyEqual(eta, this->etamax())) return 0.0;
  const double result = std::max(0.0, pow(eta, mn - 1)/(mK*this->referenceDensity()));
  ENSURE(result >= 0.0);
  return result;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
MurnahanEquationOfState<Dimension>::valid() const {
  return (SolidEquationOfState<Dimension>::valid() && 
          mn > 0.0 &&
          mK > 0.0 &&
          mAtomicWeight > 0.0 &&
          mCv > 0.0);
}

}
