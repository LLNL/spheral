//---------------------------------Spheral++----------------------------------//
// YieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent yield strength state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#include "YieldStrengthPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YieldStrengthPolicy<Dimension>::
YieldStrengthPolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
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
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::yieldStrength and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Scalar> stateFields = state.fields(fieldKey, Scalar());
  const unsigned numFields = stateFields.numFields();

  // Get the mass density, specific thermal energy, pressure,
  // plastic strain, and plastic strain rate from the state.
  const FieldList<Dimension, Scalar> massDensity = state.fields(HydroFieldNames::massDensity, 0.0);
  const FieldList<Dimension, Scalar> energy = state.fields(HydroFieldNames::specificThermalEnergy, 0.0);
  const FieldList<Dimension, Scalar> P = state.fields(HydroFieldNames::pressure, 0.0);
  const FieldList<Dimension, Scalar> PS = state.fields(SolidFieldNames::plasticStrain, 0.0);
  const FieldList<Dimension, Scalar> PSR = derivs.fields(SolidFieldNames::plasticStrainRate, 0.0);

  // Walk the individual fields.
  for (unsigned k = 0; k != numFields; ++k) {

    // Get the strength model.  This cast is ugly, but is a work-around for now.
    const SolidNodeList<Dimension>* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(stateFields[k]->nodeListPtr());
    CHECK(solidNodeListPtr != 0);
    const StrengthModel<Dimension>& strengthModel = solidNodeListPtr->strengthModel();

    // Now set the yield strength.
    strengthModel.yieldStrength(*stateFields[k], *massDensity[k], *energy[k], *P[k], *PS[k], *PSR[k]);
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

