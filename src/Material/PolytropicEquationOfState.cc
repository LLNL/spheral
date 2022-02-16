//---------------------------------Spheral++----------------------------------//
// PolytropicEquationOfState -- The gamma law gas equation of state.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//

#include "PolytropicEquationOfState.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"
#include "Utilities/SpheralFunctions.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct with the given polytropic constant, index and mean molecular 
// weight.
//------------------------------------------------------------------------------
template<typename Dimension>
PolytropicEquationOfState<Dimension>::
PolytropicEquationOfState(const double K,
                          const double index,
                          const double mu,
                          const PhysicalConstants& constants,
                          const double minimumPressure,
                          const double maximumPressure,
                          const MaterialPressureMinType minPressureType):
  EquationOfState<Dimension>(constants, minimumPressure, maximumPressure, minPressureType),
  mPolytropicConstant(K),
  mPolytropicIndex(index),
  mGamma(0.0),
  mGamma1(0.0),
  mMolecularWeight(mu),
  mExternalPressure(0.0) {
  REQUIRE(index != 0.0);
  mGamma = (index + 1.0)/index;
  mGamma1 = mGamma - 1.0;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PolytropicEquationOfState<Dimension>::
~PolytropicEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    Pressure(i) = this->pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    temperature(i) = this->temperature(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    specificThermalEnergy(i) = this->specificThermalEnergy(massDensity(i), temperature(i));
  }
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& /*massDensity*/,
                const Field<Dimension, Scalar>& /*temperature*/) const {
  CHECK(valid());
  const double kB = mConstants.kB();
  const double mp = mConstants.protonMass();
  const double Cv = kB/(mGamma1*mMolecularWeight*mp);
  specificHeat = Cv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    soundSpeed(i) = this->soundSpeed(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
              const Field<Dimension, Scalar>& /*massDensity*/,
              const Field<Dimension, Scalar>& /*specificThermalEnergy*/) const {
  CHECK(valid());
  gamma = mGamma;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a 
// polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  setPressure(bulkModulus, massDensity, specificThermalEnergy);
  bulkModulus += mExternalPressure;
  bulkModulus *= mGamma;
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PolytropicEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    entropy(i) = pressure(massDensity(i), specificThermalEnergy(i))*safeInvVar(pow(massDensity(i), mGamma));
  }
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
pressure(const Scalar massDensity,
         const Scalar /*specificThermalEnergy*/) const {
  CHECK(valid());
  return this->applyPressureLimits(mPolytropicConstant*pow(massDensity, mGamma) - mExternalPressure);
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
temperature(const Scalar /*massDensity*/,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const double kB = mConstants.kB();
  const double mp = mConstants.protonMass();
  return mGamma1*mMolecularWeight*mp/kB*specificThermalEnergy;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
specificThermalEnergy(const Scalar /*massDensity*/,
                      const Scalar temperature) const {
  CHECK(valid());
  const double kB = mConstants.kB();
  return kB/(mGamma1*mMolecularWeight)*temperature;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
specificHeat(const Scalar /*massDensity*/,
             const Scalar /*temperature*/) const {
  CHECK(valid());
  const double kB = mConstants.kB();
  const double mp = mConstants.protonMass();
  return kB/(mGamma1*mMolecularWeight*mp);
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
soundSpeed(const Scalar massDensity,
           const Scalar /*specificThermalEnergy*/) const {
  CHECK(valid());
  const double c2 = mPolytropicConstant*pow(massDensity, mGamma1);
  CHECK(c2 >= 0.0);
  return sqrt(c2);
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
gamma(const Scalar /*massDensity*/,
      const Scalar /*specificThermalEnergy*/) const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Calculate an individual bulk modulus.  
// This is just the pressure for a polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return pressure(massDensity, specificThermalEnergy) + mExternalPressure;
}

//------------------------------------------------------------------------------
// Calculate an entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension>::
entropy(const Scalar massDensity,
        const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return this->pressure(massDensity, specificThermalEnergy)*safeInvVar(pow(massDensity, mGamma));
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PolytropicEquationOfState<Dimension>::valid() const {
  return (mPolytropicConstant >= 0.0 &&
          mPolytropicIndex > 0.0 &&
          fuzzyEqual(mGamma, 1.0/mPolytropicIndex + 1.0) &&
          fuzzyEqual(mGamma1, mGamma - 1.0) &&
          mMolecularWeight > 0.0);
}

}
