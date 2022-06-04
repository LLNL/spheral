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
  mAlphaPtr(nullptr),
  mAlpha0Ptr(nullptr),
  mC0Ptr(nullptr) {
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
// Set the specific HEAT.
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
  REQUIRE(mAlpha0Ptr->nodeListPtr() == soundSpeed.nodeListPtr());
  REQUIRE(mC0Ptr->nodeListPtr() == soundSpeed.nodeListPtr());

  const auto rhoS = (*mAlphaPtr)*massDensity;
  mSolidEOS.setSoundSpeed(soundSpeed, rhoS, specificThermalEnergy);
  
  // Now apply the porosity modifier.
  const auto n = soundSpeed.numInternalElements();
  for (auto i = 0u; i != n; ++i) {
    const auto alpha0i = (*mAlpha0Ptr)(i);
    const auto alphai = (*mAlphaPtr)(i);
    const auto c0i = (*mC0Ptr)(i);
    CHECK(alpha0i >= 1.0 and alphai >= 1.0);
    CHECK(c0i > 0.0);
    soundSpeed(i) += (alphai - 1.0)*safeInv(alpha0i - 1.0)*(c0i - soundSpeed(i));
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
  
  this->setSoundSpeed(bulkModulus, massDensity, specificThermalEnergy);

  // Now apply the porosity modifier.
  const unsigned n = bulkModulus.numInternalElements();
  for (unsigned i = 0; i != n; ++i) {
    bulkModulus(i) *= massDensity(i)*bulkModulus(i);
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
// alpha
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
const Field<Dimension, typename Dimension::Scalar>&
PorousEquationOfState<Dimension>::
alpha0() const {
  return *mAlpha0Ptr;
}

template<typename Dimension>
void
PorousEquationOfState<Dimension>::
alpha0(const Field<Dimension, typename Dimension::Scalar>& x) {
  mAlpha0Ptr = &x;
}

//------------------------------------------------------------------------------
// C0
//------------------------------------------------------------------------------
template<typename Dimension>
const Field<Dimension, typename Dimension::Scalar>&
PorousEquationOfState<Dimension>::
c0() const {
  return *mC0Ptr;
}

template<typename Dimension>
void
PorousEquationOfState<Dimension>::
c0(const Field<Dimension, typename Dimension::Scalar>& x) {
  mC0Ptr = &x;
}

}
