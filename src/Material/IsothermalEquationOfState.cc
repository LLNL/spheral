//---------------------------------Spheral++----------------------------------//
// IsothermalEquationOfState
//----------------------------------------------------------------------------//
#include "IsothermalEquationOfState.hh"
#include "PhysicalConstants.hh"
#include "Field/Field.hh"
#include "Infrastructure/SpheralFunctions.hh"
#include "DBC.hh"

namespace Spheral {
namespace Material {

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
IsothermalEquationOfState<Dimension, Constants>::
IsothermalEquationOfState(const double K,
                          const double mu,
                          const double minimumPressure,
                          const double maximumPressure):
  EquationOfState<Dimension>(minimumPressure, maximumPressure),
  mK(K),
  mCs(sqrt(K)),
  mMolecularWeight(mu),
  mExternalPressure(0.0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
IsothermalEquationOfState<Dimension, Constants>::~IsothermalEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
IsothermalEquationOfState<Dimension, Constants>::
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
template<typename Dimension, typename Constants>
void
IsothermalEquationOfState<Dimension, Constants>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  temperature = 0.0;
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
IsothermalEquationOfState<Dimension, Constants>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(valid());
  specificThermalEnergy = 0.0;
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
IsothermalEquationOfState<Dimension, Constants>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(valid());
  specificHeat = 0.0;
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
IsothermalEquationOfState<Dimension, Constants>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  soundSpeed = mCs;
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
IsothermalEquationOfState<Dimension, Constants>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  gamma = 1.0;
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).  This is just the pressure for an
// isothermal gas.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
void
IsothermalEquationOfState<Dimension, Constants>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(valid());
  setPressure(bulkModulus, massDensity, specificThermalEnergy);
  bulkModulus += mExternalPressure;
}

//------------------------------------------------------------------------------
// Calculate an individual pressure.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension, Constants>::
pressure(const Scalar massDensity,
         const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return max(this->minimumPressure(), 
             min(this->maximumPressure(), 
                 mK*massDensity - mExternalPressure));
}

//------------------------------------------------------------------------------
// Calculate an individual temperature.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension, Constants>::
temperature(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return 0.0;
}

//------------------------------------------------------------------------------
// Calculate an individual specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension, Constants>::
specificThermalEnergy(const Scalar massDensity,
                      const Scalar temperature) const {
  REQUIRE(valid());
  return 0.0;
}

//------------------------------------------------------------------------------
// Calculate an individual specific heat.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension, Constants>::
specificHeat(const Scalar massDensity,
             const Scalar temperature) const {
  REQUIRE(valid());
  return 0.0;
}

//------------------------------------------------------------------------------
// Calculate an individual sound speed.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension, Constants>::
soundSpeed(const Scalar massDensity,
           const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return mCs;
}

//------------------------------------------------------------------------------
// Get gamma.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension, Constants>::
gamma(const Scalar massDensity,
      const Scalar specificThermalEnergy) const {
  return 1.0;
}

//------------------------------------------------------------------------------
// Calculate an individual bulk modulus.  
// This is just the pressure for a polytropic gas.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
typename Dimension::Scalar
IsothermalEquationOfState<Dimension, Constants>::
bulkModulus(const Scalar massDensity,
            const Scalar specificThermalEnergy) const {
  REQUIRE(valid());
  return pressure(massDensity, specificThermalEnergy) + mExternalPressure;
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension, typename Constants>
bool
IsothermalEquationOfState<Dimension, Constants>::valid() const {
  return (mK > 0.0 &&
          mCs > 0.0 &&
          mMolecularWeight > 0.0);
}

}
}

