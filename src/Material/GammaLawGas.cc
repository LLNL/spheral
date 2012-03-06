//---------------------------------Spheral++----------------------------------//
// GammaLawGas -- The gamma law gas equation of state.
//
// Created by JMO, Mon Dec  6 21:36:45 PST 1999
//----------------------------------------------------------------------------//
#include <cmath>

#include "GammaLawGas.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace Material {

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Construct with the given gamma and mu.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
GammaLawGas<Dimension, Constants>::GammaLawGas(const double gamma,
                                               const double mu,
                                               const double minimumPressure,
                                               const double maximumPressure):
  EquationOfState<Dimension>(minimumPressure, maximumPressure),
  mGamma(gamma),
  mMolecularWeight(mu) {
  mGamma1 = mGamma - 1.0;
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
GammaLawGas<Dimension, Constants>::~GammaLawGas() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
GammaLawGas<Dimension, Constants>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    Pressure(i) = pressure(massDensity(i), specificThermalEnergy(i));
  }
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
GammaLawGas<Dimension, Constants>::
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
GammaLawGas<Dimension, Constants>::
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
GammaLawGas<Dimension, Constants>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  CHECK(valid());
  double Cv = Constants::kBoltzmann/(mGamma1*mMolecularWeight*Constants::ProtonMass);
  specificHeat = Cv;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
GammaLawGas<Dimension, Constants>::
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
GammaLawGas<Dimension, Constants>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  gamma = mGamma;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for a gamma
// law gas.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
GammaLawGas<Dimension, Constants>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  setPressure(bulkModulus, massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
GammaLawGas<Dimension, Constants>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return max(this->minimumPressure(), min(this->maximumPressure(), 
                                          mGamma1*massDensity*specificThermalEnergy));
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
GammaLawGas<Dimension, Constants>::
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
GammaLawGas<Dimension, Constants>::
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
GammaLawGas<Dimension, Constants>::
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
GammaLawGas<Dimension, Constants>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return sqrt(max(0.0, mGamma*mGamma1*specificThermalEnergy));
}

//------------------------------------------------------------------------------
// Calculate an individual bulk modulus.  
// This is just the pressure for a gamma law gas.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
GammaLawGas<Dimension, Constants>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  CHECK(valid());
  return pressure(massDensity, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
GammaLawGas<Dimension, Constants>::gamma(const Scalar massDensity,
					 const Scalar specificThermalEnergy) const {
  return mGamma;
}

//------------------------------------------------------------------------------
// Get and set gamma.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
GammaLawGas<Dimension, Constants>::getGamma() const {
  return mGamma;
}

template<typename Dimension, typename Constants>
void
GammaLawGas<Dimension, Constants>::setGamma(double gamma) {
  mGamma = gamma;
  mGamma1 = mGamma - 1.0;
}

//------------------------------------------------------------------------------
// Get and set the molecular weight.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
double
GammaLawGas<Dimension, Constants>::getMolecularWeight() const {
  return mMolecularWeight;
}

template<typename Dimension, typename Constants>
void
GammaLawGas<Dimension, Constants>::setMolecularWeight(double molecularWeight) {
  mMolecularWeight = molecularWeight;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
bool
GammaLawGas<Dimension, Constants>::valid() const {
  return (mGamma > 0.0 &&
          mMolecularWeight > 0.0 &&
          mGamma1 == mGamma - 1.0);
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
    template class GammaLawGas<Dim<1>, PhysicalConstants<MKSUnits> >;
    template class GammaLawGas<Dim<2>, PhysicalConstants<MKSUnits> >;
    template class GammaLawGas<Dim<3>, PhysicalConstants<MKSUnits> >;
    template class GammaLawGas<Dim<1>, PhysicalConstants<CGSUnits> >;
    template class GammaLawGas<Dim<2>, PhysicalConstants<CGSUnits> >;
    template class GammaLawGas<Dim<3>, PhysicalConstants<CGSUnits> >;
    template class GammaLawGas<Dim<1>, PhysicalConstants<CosmologicalUnits> >;
    template class GammaLawGas<Dim<2>, PhysicalConstants<CosmologicalUnits> >;
    template class GammaLawGas<Dim<3>, PhysicalConstants<CosmologicalUnits> >;
  }
}
