//---------------------------------Spheral++----------------------------------//
// YieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent yield strength state.
//
// Created by JMO, Thu Sep 30 23:01:00 PDT 2004
//----------------------------------------------------------------------------//
#include "YieldStrengthPolicy.hh"
#include "NodeList/SolidNodeList.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Strength/SolidFieldNames.hh"
#include "Hydro/HydroFieldNames.hh"
#include "DataBase/UpdatePolicyBase.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "DataBase/IncrementState.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
YieldStrengthPolicy<Dimension>::
YieldStrengthPolicy():
  FieldUpdatePolicy<Dimension>({HydroFieldNames::massDensity,
                                HydroFieldNames::specificThermalEnergy,
                                HydroFieldNames::pressure,
                                SolidFieldNames::plasticStrain,
                                SolidFieldNames::tensorDamage,
                                IncrementState<Dimension, Scalar>::prefix() + SolidFieldNames::plasticStrain}) {
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
  REQUIRE(fieldKey == SolidFieldNames::yieldStrength);
  auto& Y = state.field(key, 0.0);

  // Get the mass density, specific thermal energy, pressure,
  // plastic strain, and plastic strain rate from the state.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& massDensity = state.field(buildKey(HydroFieldNames::massDensity), 0.0);
  const auto& energy = state.field(buildKey(HydroFieldNames::specificThermalEnergy), 0.0);
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), 0.0);
  const auto& PS = state.field(buildKey(SolidFieldNames::plasticStrain), 0.0);
  const auto& PSR = derivs.field(buildKey(SolidFieldNames::plasticStrainRate), 0.0);
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero);

  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(Y.nodeListPtr());
  CHECK(solidNodeListPtr != nullptr);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Now set the yield strength.
  strengthModel.yieldStrength(Y, massDensity, energy, P, PS, PSR, D);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
YieldStrengthPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto* rhsPtr = dynamic_cast<const YieldStrengthPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

