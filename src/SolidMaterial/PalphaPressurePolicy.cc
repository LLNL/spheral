//---------------------------------Spheral++----------------------------------//
// PalphaPressurePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the dependent pressure state for use with the P-alpha
// porosity model.
//
// Created by JMO, Wed Aug 30 13:36:37 PDT 2023
//----------------------------------------------------------------------------//

#include "PalphaPressurePolicy.hh"
#include "Hydro/HydroFieldNames.hh"
#include "FSISPH/FSIFieldNames.hh"
#include "DataBase/FieldUpdatePolicyBase.hh"
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
PalphaPressurePolicy<Dimension>::
PalphaPressurePolicy():
  FieldListUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
                                                                   HydroFieldNames::specificThermalEnergy) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PalphaPressurePolicy<Dimension>::
~PalphaPressurePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PalphaPressurePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE((fieldKey == HydroFieldNames::pressure or 
           fieldKey == FSIFieldNames::rawPressure) and
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto P = state.fields(fieldKey, Scalar());
  const auto numFields = P.numFields();

  // Get the mass density and specific thermal energy fields from the state.
  const auto massDensity = state.fields(HydroFieldNames::massDensity, Scalar());
  const auto eps = state.fields(HydroFieldNames::specificThermalEnergy, Scalar());
  CHECK(massDensity.numFields() == numFields);
  CHECK(eps.numFields() == numFields);

  // Walk the fields.
  for (auto i = 0u; i < numFields; ++i) {

    // Get the eos.  This cast is ugly, but is a work-around for now.
    const auto  fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(P[i]->nodeListPtr());
    CHECK(fluidNodeListPtr != nullptr);
    const auto& eos = fluidNodeListPtr->equationOfState();

    // Check if we need to update the pressure derivatives by seeing if they're registered for this NodeList
    const auto nodeListName = fluidNodeListPtr->name();
    const auto dPduKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialEps, nodeListName);
    const auto dPdrKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialRho, nodeListName);
    if (state.registered(dPduKey)) {
      CHECK(state.registered(dPdrKey));
      auto& dPdu = state.field(dPduKey, 0.0);
      auto& dPdr = state.field(dPdrKey, 0.0);
      eos.setPressureAndDerivs(*P[i], dPdu, dPdr, *massDensity[i], *eps[i]);
    } else {
      eos.setPressure(*P[i], *massDensity[i], *eps[i]);
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PalphaPressurePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const PalphaPressurePolicy<Dimension>*>(&rhs);
  return not (rhsPtr == nullptr);
}

}

