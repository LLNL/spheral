//---------------------------------Spheral++----------------------------------//
// PorousBulkModulusPolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent bulk modulus state in the presence of porosity.
//
// Created by JMO, Thu Dec 17 09:36:30 PST 2015
//----------------------------------------------------------------------------//

#include "PorousBulkModulusPolicy.hh"
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
PorousBulkModulusPolicy<Dimension>::
PorousBulkModulusPolicy():
  FieldUpdatePolicy<Dimension>(SolidFieldNames::porositySolidDensity,
                               HydroFieldNames::specificThermalEnergy,
                               SolidFieldNames::porosityAlpha) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PorousBulkModulusPolicy<Dimension>::
~PorousBulkModulusPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PorousBulkModulusPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == SolidFieldNames::bulkModulus);
  auto& K = state.field(key, Scalar());

  // Get the solid mass density and specific thermal energy fields from the state.
  const auto& rhoS = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porositySolidDensity, nodeListKey), Scalar());
  const auto& eps = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), Scalar());
  const auto& alpha = state.field(StateBase<Dimension>::buildFieldKey(SolidFieldNames::porosityAlpha, nodeListKey), Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(K.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Set the solid phase bulk modulus
  eos.setBulkModulus(K, rhoS, eps);

  // Scale the result by the distention
  K /= alpha;
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PorousBulkModulusPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  const auto* rhsPtr = dynamic_cast<const PorousBulkModulusPolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

