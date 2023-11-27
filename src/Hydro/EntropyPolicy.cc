//---------------------------------Spheral++----------------------------------//
// EntropyPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the entropy.
//
// Created by JMO, Sun Mar  6 21:57:47 PST 2016
//----------------------------------------------------------------------------//

#include "EntropyPolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
EntropyPolicy<Dimension>::
EntropyPolicy():
  FieldUpdatePolicy<Dimension>({HydroFieldNames::massDensity,
                                HydroFieldNames::specificThermalEnergy}) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
EntropyPolicy<Dimension>::
~EntropyPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
EntropyPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::entropy);
  auto& entropy = state.field(key, Scalar());

  // Get the mass density and specific thermal energy fields from the state.
  const auto& massDensity = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey), Scalar());
  const auto& eps = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(entropy.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Now set the entropy for this field.
  eos.setEntropy(entropy, massDensity, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
EntropyPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const EntropyPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

