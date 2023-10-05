//---------------------------------Spheral++----------------------------------//
// StrengthSoundSpeedPolicy -- Override the default sound speed policy in the 
// presence of strength.
//----------------------------------------------------------------------------//
#include "StrengthSoundSpeedPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "NodeList/SolidNodeList.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrengthSoundSpeedPolicy<Dimension>::
StrengthSoundSpeedPolicy():
  SoundSpeedPolicy<Dimension>() {
  this->addDependency(HydroFieldNames::pressure);
  this->addDependency(SolidFieldNames::tensorDamage);
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StrengthSoundSpeedPolicy<Dimension>::
~StrengthSoundSpeedPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StrengthSoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed);
  auto& cs = state.field(key, 0.0);

  // Have the base class set the initial sound speed.
  SoundSpeedPolicy<Dimension>::update(key, state, derivs, multiplier, t, dt);

  // Get the density, energy, and pressure fields from the state.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& rho = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero);

  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(cs.nodeListPtr());
  CHECK(solidNodeListPtr != 0);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Set the full sound speed.
  if (strengthModel.providesSoundSpeed()) strengthModel.soundSpeed(cs, rho, eps, P, cs, D);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StrengthSoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const StrengthSoundSpeedPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

