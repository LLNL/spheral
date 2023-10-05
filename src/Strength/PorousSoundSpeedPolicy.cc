//---------------------------------Spheral++----------------------------------//
// PorousSoundSpeedPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent sound speed state accounting for the porosity
// using a P-alpha model
//
// Created by JMO, Thu Sep 28 16:39:01 PDT 2023
//----------------------------------------------------------------------------//

#include "PorousSoundSpeedPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/SolidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousSoundSpeedPolicy<Dimension>::
PorousSoundSpeedPolicy():
  FieldUpdatePolicy<Dimension>(SolidFieldNames::porositySolidDensity,
                               HydroFieldNames::specificThermalEnergy,
                               SolidFieldNames::porosityAlpha) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousSoundSpeedPolicy<Dimension>::
~PorousSoundSpeedPolicy() {
}

//------------------------------------------------------------------------------
// Update the sound speed.
// When this is called the generic FieldList SoundSpeedPolicy will already have
// executed, so we're going to override the sound speed for our NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousSoundSpeedPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Split the field and NodeList keys.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::soundSpeed);

  // Grab the state we need.
  auto&       cs = state.field(key, 0.0);
  const auto& alpha = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha, nodeListKey), 0.0);
  const auto& alpha0 = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha0, nodeListKey), 0.0);
  const auto& cs0 = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porosityc0, nodeListKey), 0.0);
  const auto& rhoS = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porositySolidDensity, nodeListKey), 0.0);
  const auto& eps = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), 0.0);

  // Extract the solid equation of state and strength model from the NodeList.
  // This involves some ugly casting -- should revisit this design.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(cs.nodeListPtr());
  CHECK(solidNodeListPtr != 0);
  const auto& eosS = solidNodeListPtr->equationOfState();
  const auto& strengthModelS = solidNodeListPtr->strengthModel();

  // Set the solid phase sound speed
  eosS.setSoundSpeed(cs, rhoS, eps);

  // Augement with the strength model if appropriate
  if (strengthModelS.providesSoundSpeed()) {
    const auto P = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey), 0.0);
    const auto D = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::tensorDamage, nodeListKey), SymTensor::zero);
    strengthModelS.soundSpeed(cs, rhoS, eps, P, cs, D);
  }

  // Assume the sound speed varies linearly from the initial porous value to
  // the current solid value with distension.
  const auto n = cs.numInternalElements();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    const auto alpha0i = alpha0(i);
    const auto alphai = alpha(i);
    const auto cs0i = cs0(i);
    CHECK(alpha0i >= 1.0 and alphai >= 1.0);
    CHECK(cs0i > 0.0);
    cs(i) += (std::min(alphai, alpha0i) - 1.0)*safeInv(alpha0i - 1.0)*(cs0i - cs(i));
    ENSURE2(cs(i) >= 0.0, "Bad porous sound speed for " << i << " : " << cs(i));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorousSoundSpeedPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const PorousSoundSpeedPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

