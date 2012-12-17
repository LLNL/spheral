//---------------------------------Spheral++----------------------------------//
// PressurePolicy -- An implementation of UpdatePolicyBase specialized
// for the updating the node weight as a dependent quantity.
//
// Created by JMO, Tue Sep 14 22:27:08 2004
//----------------------------------------------------------------------------//

#include "PressurePolicy.hh"
#include "HydroFieldNames.hh"
#include "DataBase/FieldUpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "NodeList/FluidNodeList.hh"
#include "Material/EquationOfState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;
using FieldSpace::Field;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PressurePolicy<Dimension>::
PressurePolicy():
  FieldUpdatePolicyBase<Dimension, typename Dimension::Scalar>(HydroFieldNames::massDensity,
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
  REQUIRE(fieldKey == HydroFieldNames::pressure);

  // Get the mass density and specific thermal energy fields from the state.
  const KeyType massDensityKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::massDensity, nodeListKey);
  const KeyType energyKey = StateBase<Dimension>::buildFieldKey(HydroFieldNames::specificThermalEnergy, nodeListKey);
  CHECK(state.registered(massDensityKey));
  CHECK(state.registered(energyKey));
  Field<Dimension, Scalar>& pressure = state.field(key, Scalar());
  const Field<Dimension, Scalar>& massDensity = state.field(massDensityKey, Scalar());
  const Field<Dimension, Scalar>& energy = state.field(energyKey, Scalar());

  // Get the eos.  This cast is ugly, but is a work-around for now.
  const FluidNodeList<Dimension>* fluidNodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(pressure.nodeListPtr());
  CHECK(fluidNodeListPtr != 0);
  const Material::EquationOfState<Dimension>& eos = fluidNodeListPtr->equationOfState();

  // Now set the pressure.
  eos.setPressure(pressure, massDensity, energy);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PressurePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const PressurePolicy<Dimension>* rhsPtr = dynamic_cast<const PressurePolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

