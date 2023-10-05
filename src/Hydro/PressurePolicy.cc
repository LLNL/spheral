//---------------------------------Spheral++----------------------------------//
// PressurePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the node weight as a dependent quantity.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//

#include "PressurePolicy.hh"
#include "HydroFieldNames.hh"
#include "FSISPH/FSIFieldNames.hh"
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
PressurePolicy<Dimension>::
PressurePolicy():
  FieldUpdatePolicy<Dimension>(HydroFieldNames::massDensity,
                               HydroFieldNames::specificThermalEnergy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PressurePolicy<Dimension>::
~PressurePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PressurePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE((fieldKey == HydroFieldNames::pressure or 
           fieldKey == FSIFieldNames::rawPressure));
  auto& P = state.field(key, Scalar());

  // Get the mass density and specific thermal energy fields from the state.
  const auto& massDensity = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey), Scalar());
  const auto& eps = state.field(State<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(P.nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Now set the pressure for this field.
  eos.setPressure(P, massDensity, eps);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PressurePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const PressurePolicy<Dimension>*>(&rhs);
  return (rhsPtr != nullptr);
}

}

