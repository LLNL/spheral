//---------------------------------Spheral++----------------------------------//
// PolytropicEquationOfState -- The gamma law gas equation of state.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//

#include "PolytropicEquationOfState.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "DBC.hh"

namespace Spheral {
namespace Material {

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Construct with the given polytropic constant, index and mean molecular 
// weight.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
PolytropicEquationOfState<Dimension, Constants>::
PolytropicEquationOfState(const double K,
                          const double index,
                          const double mu,
                          const double minimumPressure,
                          const double maximumPressure):
  EquationOfState<Dimension>(minimumPressure, maximumPressure),
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
template<typename Dimension, typename Constants>
PolytropicEquationOfState<Dimension, Constants>::~PolytropicEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
PolytropicEquationOfState<Dimension, Constants>::
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
template<typename Dimension, typename Constants>
void
PolytropicEquationOfState<Dimension, Constants>::
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
template<typename Dimension, typename Constants>
void
PolytropicEquationOfState<Dimension, Constants>::
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
template<typename Dimension, typename Constants>
void
PolytropicEquationOfState<Dimension, Constants>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
  const double Cv = Constants::kBoltzmann/(mGamma1*mMolecularWeight*Constants::ProtonMass);
  specificHeat = Cv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
PolytropicEquationOfState<Dimension, Constants>::
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
template<typename Dimension, typename Constants>
void
PolytropicEquationOfState<Dimension, Constants>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  gamma = mGamma;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a 
// polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
PolytropicEquationOfState<Dimension, Constants>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  setPressure(bulkModulus, massDensity, specificThermalEnergy);
  bulkModulus += mExternalPressure;
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension, Constants>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return max(this->minimumPressure(), 
             min(this->maximumPressure(), 
                 mPolytropicConstant*pow(massDensity, mGamma) - mExternalPressure));
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension, Constants>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return mGamma1*mMolecularWeight*Constants::ProtonMass/Constants::kBoltzmann*specificThermalEnergy;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension, Constants>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  CHECK(valid());
  return Constants::kBoltzmann/(mGamma1*mMolecularWeight)*temperature;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension, Constants>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  CHECK(valid());
  return Constants::kBoltzmann/(mGamma1*mMolecularWeight*Constants::ProtonMass);
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension, Constants>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  CHECK(valid());
  const double c2 = mPolytropicConstant*pow(massDensity, mGamma1);
  CHECK(c2 >= 0.0);
  return sqrt(c2);
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension, Constants>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Calculate an individual bulk modulus.  
// This is just the pressure for a polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
PolytropicEquationOfState<Dimension, Constants>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return pressure(massDensity, specificThermalEnergy) + mExternalPressure;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
bool
PolytropicEquationOfState<Dimension, Constants>::valid() const {
  return (mPolytropicConstant >= 0.0 &&
          mPolytropicIndex > 0.0 &&
          fuzzyEqual(mGamma, 1.0/mPolytropicIndex + 1.0) &&
          fuzzyEqual(mGamma1, mGamma - 1.0) &&
          mMolecularWeight > 0.0);
}

}
}

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PhysicalConstants.hh"
#include "MKSUnits.hh"
#include "CGSUnits.hh"
#include "CosmologicalUnits.hh"
namespace Spheral {
  namespace Material {
    template class PolytropicEquationOfState<Dim<1>, PhysicalConstants<MKSUnits> >;
    template class PolytropicEquationOfState<Dim<2>, PhysicalConstants<MKSUnits> >;
    template class PolytropicEquationOfState<Dim<3>, PhysicalConstants<MKSUnits> >;
    template class PolytropicEquationOfState<Dim<1>, PhysicalConstants<CGSUnits> >;
    template class PolytropicEquationOfState<Dim<2>, PhysicalConstants<CGSUnits> >;
    template class PolytropicEquationOfState<Dim<3>, PhysicalConstants<CGSUnits> >;
    template class PolytropicEquationOfState<Dim<1>, PhysicalConstants<CosmologicalUnits> >;
    template class PolytropicEquationOfState<Dim<2>, PhysicalConstants<CosmologicalUnits> >;
    template class PolytropicEquationOfState<Dim<3>, PhysicalConstants<CosmologicalUnits> >;
  }
}
