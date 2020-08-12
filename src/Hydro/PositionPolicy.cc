//---------------------------------Spheral++----------------------------------//
// PositionPolicy -- An implementation of UpdatePolicyBase specialized for the
// updating the position.
//
// This version ignores the XSPH approximation in order to time center the
// velocity for updating the position.  This is intended for use with the 
// compatible energy evolution hydro approximation.
//
// Created by JMO, Mon Jun 19 22:06:07 PDT 2006
//----------------------------------------------------------------------------//
#include "PositionPolicy.hh"
#include "HydroFieldNames.hh"
#include "NodeList/FluidNodeList.hh"
#include "DataBase/FieldListUpdatePolicyBase.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Geometry/Dimension.hh"

namespace Spheral {


//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PositionPolicy<Dimension>::
PositionPolicy():
  IncrementFieldList<Dimension, typename Dimension::Vector>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PositionPolicy<Dimension>::
~PositionPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
PositionPolicy<Dimension>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& /*derivs*/,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::position and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  auto r = state.fields(fieldKey, Vector::zero);
  const auto numFields = r.numFields();

  // Get the velocity and acceleration fields.
  const auto vel = state.fields(HydroFieldNames::velocity, Vector::zero);
  // const auto dvel = derivs.fields(HydroFieldNames::hydroAcceleration, Vector::zero);

  // Walk the fields.
  for (auto i = 0u; i != numFields; ++i) {
    const auto n = r[i]->numInternalElements();
    for (auto j = 0u; j < n; ++j) {
      r(i,j) += multiplier*vel(i,j);
      // r(i,j) += multiplier*(vel(i,j) + 0.5*multiplier*dvel(i,j));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
PositionPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const PositionPolicy<Dimension>* rhsPtr = dynamic_cast<const PositionPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

