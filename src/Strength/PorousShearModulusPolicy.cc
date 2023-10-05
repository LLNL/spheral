//---------------------------------Spheral++----------------------------------//
// PorousShearModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent gamma state in the presence of porosity.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//

#include "PorousShearModulusPolicy.hh"
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
PorousShearModulusPolicy<Dimension>::
PorousShearModulusPolicy():
  FieldUpdatePolicy<Dimension>(SolidFieldNames::porositySolidDensity,
                               HydroFieldNames::specificThermalEnergy,
                               HydroFieldNames::pressure,
                               SolidFieldNames::tensorDamage,
                               SolidFieldNames::porosityAlpha) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousShearModulusPolicy<Dimension>::
~PorousShearModulusPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousShearModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::shearModulus);
  auto& mu = state.field(key, Scalar());

  // Get the strength model.  This cast is ugly, but is a work-around for now.
  const auto* solidNodeListPtr = dynamic_cast<const SolidNodeList<Dimension>*>(mu->nodeListPtr());
  CHECK(solidNodeListPtr != nullptr);
  const auto& strengthModel = solidNodeListPtr->strengthModel();

  // Get the state fields we need
  const auto& rhoS = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porositySolidDensity, nodeListKey), Scalar());
  const auto& eps = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), Scalar());
  const auto& P = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey), Scalar());
  const auto& D = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::tensorDamage, nodeListKey), SymTensor::zero);
  const auto& alpha = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha, nodeListKey), Scalar());

  // We actually need the solid phase pressure
  const auto PS = P*alpha;

  // Set the solid phase shear modulus
  strengthModel.shearModulus(mu, rhoS, eps, PS, D);

  // Scale the result by the distention
  mu /= alpha;
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorousShearModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  const auto* rhsPtr = dynamic_cast<const PorousShearModulusPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

