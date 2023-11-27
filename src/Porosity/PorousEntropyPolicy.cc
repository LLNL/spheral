//---------------------------------Spheral++----------------------------------//
// PorousEntropyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent entropy state accounting for the porosity
// using a P-alpha model
//
// Created by JMO, Thu Oct  5 11:29:39 PDT 2023
//----------------------------------------------------------------------------//

#include "PorousEntropyPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousEntropyPolicy<Dimension>::
PorousEntropyPolicy():
  FieldUpdatePolicy<Dimension>({SolidFieldNames::porositySolidDensity,
                                HydroFieldNames::specificThermalEnergy,
                                SolidFieldNames::porosityAlpha}) {
}

//------------------------------------------------------------------------------
// Update the sound speed.
// When this is called the generic FieldList EntropyPolicy will already have
// executed, so we're going to override the sound speed for our NodeList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousEntropyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Split the field and NodeList keys.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::entropy);

  // Grab the state we need.
  auto&       entropy = state.field(key, 0.0);
  const auto& rhoS = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porositySolidDensity, nodeListKey), 0.0);
  const auto& eps = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), 0.0);

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(entropy.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Set the solid phase entropy
  eos.setEntropy(entropy, rhoS, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorousEntropyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const PorousEntropyPolicy<Dimension>*>(&rhs) != nullptr;
}

}

