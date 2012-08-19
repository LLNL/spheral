//---------------------------------Spheral++----------------------------------//
// StrainPorosity -- Strain porosity EOS modifier.
// 
// See header for references and such.
//----------------------------------------------------------------------------//

#include "PorousEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {
namespace SolidMaterial {

using namespace std;
using std::min;
using std::max;
using std::abs;

using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousEquationOfState<Dimension>::
PorousEquationOfState(const SolidEquationOfState<Dimension>& solidEOS):
  SolidEquationOfState<Dimension>(solidEOS.referenceDensity(),
                                  solidEOS.etamin(),
                                  solidEOS.etamax(),
                                  solidEOS.minimumPressure(),
                                  solidEOS.maximumPressure()),
  mSolidEOS(solidEOS),
  mAlphaPtr(0) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousEquationOfState<Dimension>::
~PorousEquationOfState() {
}

//------------------------------------------------------------------------------
// Set the pressure.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setPressure(Field<Dimension, Scalar>& Pressure,
            const Field<Dimension, Scalar>& massDensity,
            const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(this->valid());
  REQUIRE(massDensity.nodeListPtr() == Pressure.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == Pressure.nodeListPtr());
  REQUIRE(mAlphaPtr->nodeListPtr() == Pressure.nodeListPtr());

  // The base EOS set's the solid (compacted) pressure.
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setPressure(Pressure, rhoS, specificThermalEnergy);

  // Now apply the porosity modifier.
  Pressure /= *mAlphaPtr;
}

//------------------------------------------------------------------------------
// Set the temperature.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setTemperature(Field<Dimension, Scalar>& temperature,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(this->valid());
  REQUIRE(massDensity.nodeListPtr() == temperature.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == temperature.nodeListPtr());
  REQUIRE(mAlphaPtr->nodeListPtr() == temperature.nodeListPtr());
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setTemperature(temperature, rhoS, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Set the specific thermal energy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setSpecificThermalEnergy(Field<Dimension, Scalar>& specificThermalEnergy,
                         const Field<Dimension, Scalar>& massDensity,
                         const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(this->valid());
  REQUIRE(massDensity.nodeListPtr() == specificThermalEnergy.nodeListPtr());
  REQUIRE(temperature.nodeListPtr() == specificThermalEnergy.nodeListPtr());
  REQUIRE(mAlphaPtr->nodeListPtr() == specificThermalEnergy.nodeListPtr());
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setSpecificThermalEnergy(specificThermalEnergy, rhoS, temperature);
}

//------------------------------------------------------------------------------
// Set the specific heat.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setSpecificHeat(Field<Dimension, Scalar>& specificHeat,
                const Field<Dimension, Scalar>& massDensity,
                const Field<Dimension, Scalar>& temperature) const {
  REQUIRE(this->valid());
  REQUIRE(massDensity.nodeListPtr() == specificHeat.nodeListPtr());
  REQUIRE(temperature.nodeListPtr() == specificHeat.nodeListPtr());
  REQUIRE(mAlphaPtr->nodeListPtr() == specificHeat.nodeListPtr());
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setSpecificHeat(specificHeat, rhoS, temperature);
}

//------------------------------------------------------------------------------
// Set the sound speed.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setSoundSpeed(Field<Dimension, Scalar>& soundSpeed,
              const Field<Dimension, Scalar>& massDensity,
              const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(this->valid());
  REQUIRE(massDensity.nodeListPtr() == soundSpeed.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == soundSpeed.nodeListPtr());
  REQUIRE(mAlphaPtr->nodeListPtr() == soundSpeed.nodeListPtr());
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setSoundSpeed(soundSpeed, rhoS, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Set gamma (ratio of specific heats).
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setGammaField(Field<Dimension, Scalar>& gamma,
	      const Field<Dimension, Scalar>& massDensity,
	      const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(this->valid());
  REQUIRE(massDensity.nodeListPtr() == gamma.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == gamma.nodeListPtr());
  REQUIRE(mAlphaPtr->nodeListPtr() == gamma.nodeListPtr());
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setGammaField(gamma, rhoS, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Set the bulk modulus (rho DP/Drho).
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setBulkModulus(Field<Dimension, Scalar>& bulkModulus,
               const Field<Dimension, Scalar>& massDensity,
               const Field<Dimension, Scalar>& specificThermalEnergy) const {
  REQUIRE(this->valid());
  REQUIRE(massDensity.nodeListPtr() == bulkModulus.nodeListPtr());
  REQUIRE(specificThermalEnergy.nodeListPtr() == bulkModulus.nodeListPtr());
  REQUIRE(mAlphaPtr->nodeListPtr() == bulkModulus.nodeListPtr());
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setBulkModulus(bulkModulus, rhoS, specificThermalEnergy);
}

//------------------------------------------------------------------------------
// Determine if the EOS is in a valid state.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorousEquationOfState<Dimension>::valid() const {
  return (mSolidEOS.valid() and mAlphaPtr != 0);
}

//------------------------------------------------------------------------------
// Access the underlying solid EOS.
//------------------------------------------------------------------------------
template<typename Dimension>
const SolidEquationOfState<Dimension>&
PorousEquationOfState<Dimension>::
solidEOS() const {
  return mSolidEOS;
}

//------------------------------------------------------------------------------
// Access the alpha field.
//------------------------------------------------------------------------------
template<typename Dimension>
const FieldSpace::Field<Dimension, typename Dimension::Scalar>&
PorousEquationOfState<Dimension>::
alpha() const {
  return *mAlphaPtr;
}

template<typename Dimension>
void
PorousEquationOfState<Dimension>::
alpha(const FieldSpace::Field<Dimension, typename Dimension::Scalar>& x) {
  mAlphaPtr = &x;
}

}
}

