//---------------------------------Spheral++----------------------------------//
// YieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent yield strength state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#include "YieldStrengthPolicy.hh"
#include "SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using FieldSpace::Field;
using NodeSpace::NodeList;
using SolidMaterial::SolidNodeList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YieldStrengthPolicy<Dimension>::
YieldStrengthPolicy():
UpdatePolicyBase<Dimension>(HydroFieldNames::massDensity,
                            HydroFieldNames::specificThermalEnergy,
                            HydroFieldNames::pressure,
                            SolidFieldNames::plasticStrain,
                            IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::plasticStrain) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YieldStrengthPolicy<Dimension>::
~YieldStrengthPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
YieldStrengthPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::yieldStrength);
  Field<Dimension, Scalar>& stateField = state.field(key, 0.0);

  // Get the mass density, specific thermal energy, pressure,
  // plastic strain, and plastic strain rate from the state.
  const KeyType massDensityKey = State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
  const KeyType energyKey = State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  const KeyType PKey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType PSKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrain, nodeListKey);
  const KeyType PSRKey = State<Dimension>::buildFieldKey(SolidFieldNames::plasticStrainRate, nodeListKey);
  CHECK(state.registered(massDensityKey));
  CHECK(state.registered(energyKey));
  CHECK(state.registered(PKey));
  CHECK(state.registered(PSKey));
  CHECK(derivs.registered(PSRKey));
  const Field<Dimension, Scalar>& massDensity = state.field(massDensityKey, 0.0);
  const Field<Dimension, Scalar>& energy = state.field(energyKey, 0.0);
  const Field<Dimension, Scalar>& P = state.field(PKey, 0.0);
  const Field<Dimension, Scalar>& PS = state.field(PSKey, 0.0);
  const Field<Dimension, Scalar>& PSR = derivs.field(PSRKey, 0.0);

  // We only do this if this is a solid node list.
  const SolidNodeList<Dimension>* nodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(stateField.nodeListPtr());
  CHECK(nodeListPtr != 0);

  // Get the strength model.
  const SolidMaterial::StrengthModel<Dimension>& strengthModel = nodeListPtr->strengthModel();

  // Now set the yield strength.
  for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {
    stateField(i) = strengthModel.yieldStrength(massDensity(i), energy(i), P(i),
                                                PS(i), PSR(i));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
YieldStrengthPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const YieldStrengthPolicy<Dimension>* rhsPtr = dynamic_cast<const YieldStrengthPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

