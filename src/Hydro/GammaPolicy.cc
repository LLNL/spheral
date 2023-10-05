//---------------------------------Spheral++----------------------------------//
// GammaPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent gamma state.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//

#include "GammaPolicy.hh"
#include "HydroFieldNames.hh"
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
GammaPolicy<Dimension>::
GammaPolicy():
  FieldUpdatePolicy<Dimension>(HydroFieldNames::massDensity,
                               HydroFieldNames::specificThermalEnergy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
GammaPolicy<Dimension>::
~GammaPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
GammaPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::gamma);
  auto& gamma = state.field(key, Scalar());

  // Get the mass density and specific thermal energy fields from the state.
  const auto& massDensity = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey), Scalar());
  const auto& eps = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(gamma.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Now set the gamma for this field.
  eos.setGammaField(gamma, massDensity, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
GammaPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const GammaPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

