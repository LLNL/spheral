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
  FieldUpdatePolicy<Dimension>(SolidFieldNames::porsitySolidDensity,
                               HydroFieldNames::specificThermalEnergy,
                               SolidFieldNames::porosityAlpha) {
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
           fieldKey == FSIFieldNames::rawPressure));
  auto& P = state.field(key, Scalar());

  // Get the solid mass density, specific thermal energy, and distention fields from the state.
  const auto& rhoS = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey), Scalar());
  const auto& eps = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey), Scalar());
  const auto& alpha = state.field(StateBase<Dimension>::buildFieldKey(HydroFieldNames::porosityAlpha, nodeListKey), Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const auto* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(P->nodeListPtr());
  CHECK(fluidNodeListPtr != nullptr);
  const auto& eos = fluidNodeListPtr->equationOfState();

  // Check if we need to update the pressure derivatives by seeing if they're registered for this NodeList
  const auto dPduKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialEps, nodeListKey);
  const auto dPdrKey = State<Dimension>::buildFieldKey(HydroFieldNames::partialPpartialRho, nodeListKey);
  if (state.registered(dPduKey)) {
    CHECK(state.registered(dPdrKey));
    auto& dPdu = state.field(dPduKey, 0.0);
    auto& dPdr = state.field(dPdrKey, 0.0);
    eos.setPressureAndDerivs(P, dPdu, dPdr, rhoS, eps);
  } else {
    eos.setPressure(P, rhoS, eps);
  }

  // Scale the pressure to the bulk value.
  P /= alpha;
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

