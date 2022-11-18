//---------------------------------Spheral++----------------------------------//
// IsothermalEquationOfState
//----------------------------------------------------------------------------//
#include "IsothermalEquationOfState.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
IsothermalEquationOfState<Dimension>::
IsothermalEquationOfState(const double K,
                          const double mu,
                          const PhysicalConstants& constants,
                          const double minimumPressure,
                          const double maximumPressure,
                          const MaterialPressureMinType minPressureType,
                          const double externalPressure):
  EquationOfState<Dimension>(constants, minimumPressure, maximumPressure, minPressureType, externalPressure),
  mK(K),
  mCs(sqrt(K)),
  mMolecularWeight(mu) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
IsothermalEquationOfState<Dimension>::
~IsothermalEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  REQUIRE(Pressure.nodeListPtr() == massDensity.nodeListPtr());
  REQUIRE(Pressure.nodeListPtr() == specificThermalEnergy.nodeListPtr());
  const int n = Pressure.nodeListPtr()->numInternalNodes();
  for (int i = 0; i != n; ++i) {
    Pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& /*massDensity*/,
               const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  temperature = 0.0;
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& /*massDensity*/,
                         const Field<Dimension, Scalar>& /*temperature*/) const {
  REQUIRE(valid());
  specificThermalEnergy = 0.0;
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  REQUIRE(valid());
  specificHeat = 0.0;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& /*massDensity*/,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  soundSpeed = mCs;
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
              const Field<Dimension, Scalar>& /*massDensity*/,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  gamma = 1.0;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for an
// isothermal gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    bulkModulus(i) = this->bulkModulus(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
IsothermalEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    entropy(i) = pressure(massDensity(i), specificThermalEnergy(i))*safeInvVar(massDensity(i));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
pressure(const Scalar massDensity,
         const Scalar /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  return this->applyPressureLimits(mK*massDensity);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
temperature(const Scalar /*massDensity*/,
            const Scalar /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  return 0.0;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
specificThermalEnergy(const Scalar /*massDensity*/,
                      const Scalar /*temperature*/) const {
  REQUIRE(valid());
  return 0.0;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
specificHeat(const Scalar /*massDensity*/,
             const Scalar /*temperature*/) const {
  REQUIRE(valid());
  return 0.0;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
soundSpeed(const Scalar /*massDensity*/,
           const Scalar /*specificThermalEnergy*/) const {
  REQUIRE(valid());
  return mCs;
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
gamma(const Scalar /*massDensity*/,
      const Scalar /*specificThermalEnergy*/) const {
  return 1.0;
}

//------------------------------------------------------------------------------
// Calculate an individual bulk modulus.  
// This is just the pressure for a polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return pressure(massDensity, specificThermalEnergy) + this->externalPressure();
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return this->pressure(massDensity, specificThermalEnergy)*safeInvVar(massDensity);
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
IsothermalEquationOfState<Dimension>::valid() const {
  return (mK > 0.0 &&
          mCs > 0.0 &&
          mMolecularWeight > 0.0);
}

}
