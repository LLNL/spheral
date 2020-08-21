//---------------------------------Spheral++----------------------------------//
// SpecificThermalEnergyVolumePolicy
//
// Created by JMO, Sun Aug 18 19:21:27 PDT 2013
//----------------------------------------------------------------------------//

#include "SpecificThermalEnergyVolumePolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Hydro/HydroFieldNames.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension>
SpecificThermalEnergyVolumePolicy<Dimension>::
SpecificThermalEnergyVolumePolicy():
  FieldUpdatePolicyBase<Dimension, Scalar>(HydroFieldNames::volume) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
SpecificThermalEnergyVolumePolicy<Dimension>::
~SpecificThermalEnergyVolumePolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
SpecificThermalEnergyVolumePolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::specificThermalEnergy);

  // Get the state fields.
  const KeyType volKey = State<Dimension>::buildFieldKey(HydroFieldNames::volume, nodeListKey);
  const KeyType vol0Key = State<Dimension>::buildFieldKey(HydroFieldNames::volume + "0", nodeListKey);
  const KeyType Pkey = State<Dimension>::buildFieldKey(HydroFieldNames::pressure, nodeListKey);
  const KeyType Qkey = State<Dimension>::buildFieldKey(HydroFieldNames::maxViscousPressure, nodeListKey);
  CHECK(state.registered(volKey));
  CHECK(state.registered(vol0Key));
  CHECK(state.registered(Pkey));
  CHECK(derivs.registered(Qkey));
  Field<Dimension, Scalar>& eps = state.field(key, 0.0);
  const Field<Dimension, Scalar>& vol = state.field(volKey, 0.0);
  const Field<Dimension, Scalar>& vol0 = state.field(vol0Key, 0.0);
  const Field<Dimension, Scalar>& P = state.field(Pkey, 0.0);
  const Field<Dimension, Scalar>& Q = derivs.field(Qkey, 0.0);

  // Do the deed.
  for (unsigned i = 0; i != eps.numInternalElements(); ++i) {
    CHECK(vol0(i) > 0.0);
    eps(i) += multiplier*((P(i) + Q(i))*(vol0(i) - vol(i)));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
SpecificThermalEnergyVolumePolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const SpecificThermalEnergyVolumePolicy<Dimension>* rhsPtr = dynamic_cast<const SpecificThermalEnergyVolumePolicy<Dimension>*>(&rhs);
  return rhsPtr != 0;
}

}

