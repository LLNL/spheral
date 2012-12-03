//---------------------------------Spheral++----------------------------------//
// StrengthModel -- The interface base class for strength models.
//
// Created by JMO, Wed Sep 8 15:18:44 2004
//----------------------------------------------------------------------------//

#include "StrengthModel.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace SolidMaterial {

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Compute an individual full sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
StrengthModel<Dimension>::
soundSpeed(const double density,
           const double specificThermalEnergy,
           const double pressure,
           const double fluidSoundSpeed) const {
  REQUIRE(distinctlyGreaterThan(density, 0.0));
  const double cs2 = fluidSoundSpeed*fluidSoundSpeed + 4.0/3.0 * shearModulus(density, specificThermalEnergy, pressure) / density;
  ENSURE(cs2 > 0.0);
  return std::sqrt(cs2);
}

//------------------------------------------------------------------------------
// Shear modulus -- field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
StrengthModel<Dimension>::
shearModulus(Field<Dimension, typename Dimension::Scalar>& shearModulus,
             const Field<Dimension, typename Dimension::Scalar>& density,
             const Field<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
             const Field<Dimension, typename Dimension::Scalar>& pressure) const {
  REQUIRE(density.nodeListPtr() == shearModulus.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == shearModulus.nodeListPtr());
  REQUIRE(pressure.nodeListPtr() == shearModulus.nodeListPtr());
  const unsigned n = shearModulus.size();
  for (unsigned i = 0; i != n; ++i) {
    shearModulus(i) = this->shearModulus(density(i), specificThermalEnergy(i), pressure(i));
  }
}

//------------------------------------------------------------------------------
// Yield strength -- field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
StrengthModel<Dimension>::
yieldStrength(Field<Dimension, typename Dimension::Scalar>& yieldStrength,
              const Field<Dimension, typename Dimension::Scalar>& density,
              const Field<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
              const Field<Dimension, typename Dimension::Scalar>& pressure,
              const Field<Dimension, typename Dimension::Scalar>& plasticStrain,
              const Field<Dimension, typename Dimension::Scalar>& plasticStrainRate) const {
  REQUIRE(density.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(pressure.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(plasticStrain.nodeListPtr() == yieldStrength.nodeListPtr());
  REQUIRE(plasticStrainRate.nodeListPtr() == yieldStrength.nodeListPtr());
  const unsigned n = yieldStrength.size();
  for (unsigned i = 0; i != n; ++i) {
    yieldStrength(i) = this->yieldStrength(density(i), specificThermalEnergy(i), pressure(i), plasticStrain(i), plasticStrainRate(i));
  }
}

//------------------------------------------------------------------------------
// Sound speed -- field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
StrengthModel<Dimension>::
soundSpeed(Field<Dimension, typename Dimension::Scalar>& soundSpeed,
           const Field<Dimension, typename Dimension::Scalar>& density,
           const Field<Dimension, typename Dimension::Scalar>& specificThermalEnergy,
           const Field<Dimension, typename Dimension::Scalar>& pressure,
           const Field<Dimension, typename Dimension::Scalar>& fluidSoundSpeed) const {
  REQUIRE(density.nodeListPtr() == soundSpeed.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == soundSpeed.nodeListPtr());
  REQUIRE(pressure.nodeListPtr() == soundSpeed.nodeListPtr());
  REQUIRE(fluidSoundSpeed.nodeListPtr() == soundSpeed.nodeListPtr());
  const unsigned n = soundSpeed.size();
  for (unsigned i = 0; i != n; ++i) {
    soundSpeed(i) = this->soundSpeed(density(i), specificThermalEnergy(i), pressure(i), fluidSoundSpeed(i));
  }
}

}
}
