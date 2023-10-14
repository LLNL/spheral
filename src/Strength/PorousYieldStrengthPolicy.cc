//---------------------------------Spheral++----------------------------------//
// PorousYieldStrengthPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent shear modulus state in the presence of porosity.
//
// Created by JMO, Thu Oct  5 11:04:29 PDT 2023
//----------------------------------------------------------------------------//

#include "PorousYieldStrengthPolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Strength/SolidFieldNames.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "NodeList/SolidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "SolidMaterial/StrengthModel.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousYieldStrengthPolicy<Dimension>::
PorousYieldStrengthPolicy():
  FieldUpdatePolicy<Dimension>(SolidFieldNames::porositySolidDensity,
                               HydroFieldNames::specificThermalEnergy,
                               SolidFieldNames::porosityAlpha) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousYieldStrengthPolicy<Dimension>::
~PorousYieldStrengthPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousYieldStrengthPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::yieldStrength);
  auto& Y = state.field(key, Scalar());

  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(Y.nodeListPtr());
  CHECK(solidNodeListPtr != nullptr);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Get the state we need
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& rhoS = state.field(buildKey(SolidFieldNames::porositySolidDensity), Scalar());
  const auto& eps = state.field(buildKey(HydroFieldNames::specificThermalEnergy), Scalar());
  const auto& P = state.field(buildKey(HydroFieldNames::pressure), Scalar());
  const auto& plasticStrain = state.field(buildKey(SolidFieldNames::plasticStrain), Scalar());
  const auto& plasticStrainRate = derivs.field(buildKey(SolidFieldNames::plasticStrainRate), Scalar());
  const auto& D = state.field(buildKey(SolidFieldNames::tensorDamage), SymTensor::zero);
  const auto& alpha = state.field(buildKey(SolidFieldNames::porosityAlpha), Scalar());

  // We actually need the solid phase pressure
  const auto PS = P*alpha;

  // Set the solid phase yield strength
  strengthModel.yieldStrength(Y, rhoS, eps, PS, plasticStrain, plasticStrainRate, D);

  // Scale the result by the distention
  Y /= alpha;
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorousYieldStrengthPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  const auto* rhsPtr = dynamic_cast<const PorousYieldStrengthPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

