//---------------------------------Spheral++----------------------------------//
// StrainPorosity -- Strain porosity EOS modifier.
// 
// See header for references and such.
//----------------------------------------------------------------------------//
#include "PorousEquationOfState.hh"
#include "Field/Field.hh"

namespace Spheral {

using std::min;
using std::max;
using std::abs;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousEquationOfState<Dimension>::
PorousEquationOfState(const EquationOfState<Dimension>& solidEOS):
  EquationOfState<Dimension>(solidEOS.constants(),
                                       solidEOS.minimumPressure(),
                                       solidEOS.maximumPressure(),
                                       solidEOS.minimumPressureType()),
  mSolidEOS(solidEOS),
  mAlphaPtr(0),
  mAlpha0(0.0),
  mC0(0.0) {
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
  const unsigned n = Pressure.numInternalElements();
  for (unsigned i = 0; i != n; ++i) {
    Pressure(i) = max(this->minimumPressure(),
                      min(this->maximumPressure(),
                          Pressure(i)/(*mAlphaPtr)(i)));
  }
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
// Based on the relation in 
// Collins, G. S., Melosh, H. J., & Wünnemann, K. (2011). Improvements to the ɛ-α porous compaction model for simulating impacts into
// high-porosity solar system objects. International Journal of Impact Engineering, 38(6), 434–439. 
// http://doi.org/10.1016/j.ijimpeng.2010.10.013
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
  REQUIRE2(mAlpha0 >= 1.0 and mC0 > 0.0, mAlpha0 << " " << mC0);
  const Field<Dimension, Scalar> rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setSoundSpeed(soundSpeed, rhoS, specificThermalEnergy);
  
  // Now apply the porosity modifier.
  const unsigned n = soundSpeed.numInternalElements();
  for (unsigned i = 0; i != n; ++i) {
    soundSpeed(i) += ((*mAlphaPtr)(i) - 1.0)*safeInv(mAlpha0 - 1.0)*(mC0 - soundSpeed(i));
  }
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

  // Now apply the porosity modifier.
  const unsigned n = bulkModulus.numInternalElements();
  for (unsigned i = 0; i != n; ++i) {
    bulkModulus(i) /= (*mAlphaPtr)(i);
  }
}

//------------------------------------------------------------------------------
// Set the entropy.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEquationOfState<Dimension>::
setEntropy(Field<Dimension, Scalar>& entropy,
           const Field<Dimension, Scalar>& massDensity,
           const Field<Dimension, Scalar>& specificThermalEnergy) const {
  CHECK(valid());
  Field<Dimension, Scalar> gamma("gamma", entropy.nodeList());
  this->setPressure(entropy, massDensity, specificThermalEnergy);
  this->setGammaField(gamma, massDensity, specificThermalEnergy);
  for (size_t i = 0; i != massDensity.numElements(); ++i) {
    entropy(i) *= safeInvVar(pow(massDensity(i), gamma(i)));
  }
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
const EquationOfState<Dimension>&
PorousEquationOfState<Dimension>::
solidEOS() const {
  return mSolidEOS;
}

//------------------------------------------------------------------------------
// Access the alpha field.
//------------------------------------------------------------------------------
template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
PorousEquationOfState<Dimension>::
alpha() const {
  return *mAlphaPtr;
}

template<typename Dimension>
void
PorousEquationOfState<Dimension>::
alpha(const Field<Dimension, typename Dimension::Scalar>& x) {
  mAlphaPtr = &x;
}

//------------------------------------------------------------------------------
// alpha0
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PorousEquationOfState<Dimension>::
alpha0() const {
  return mAlpha0;
}

template<typename Dimension>
void
PorousEquationOfState<Dimension>::
alpha0(typename Dimension::Scalar x) {
  mAlpha0 = x;
}

//------------------------------------------------------------------------------
// c0
//------------------------------------------------------------------------------
template<typename Dimension>
typename Dimension::Scalar
PorousEquationOfState<Dimension>::
c0() const {
  return mC0;
}

template<typename Dimension>
void
PorousEquationOfState<Dimension>::
c0(typename Dimension::Scalar x) {
  mC0 = x;
}

}
