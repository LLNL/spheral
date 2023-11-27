//---------------------------------Spheral++----------------------------------//
// PorousGammaPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent gamma state in the presence of porosity.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//

#include "PorousGammaPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
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
PorousGammaPolicy<Dimension>::
PorousGammaPolicy():
  FieldUpdatePolicy<Dimension>({SolidFieldNames::porositySolidDensity,
                                HydroFieldNames::specificThermalEnergy,
                                SolidFieldNames::porosityAlpha}) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousGammaPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::gamma);
  auto& gamma = state.field(key, Scalar());

  // Get the solid mass density and specific thermal energy fields from the state.
  const auto& rhoS = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porositySolidDensity, nodeListKey), Scalar());
  const auto& eps = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(gamma.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Now set the gamma for this field.
  eos.setGammaField(gamma, rhoS, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorousGammaPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const PorousGammaPolicy<Dimension>*>(&rhs) != nullptr;
}

}

