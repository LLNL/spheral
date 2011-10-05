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
#include "DataBase/FieldUpdatePolicyBase.hh"
#include "DataBase/IncrementState.hh"
#include "DataBase/ReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using namespace std;
using FieldSpace::Field;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;

//------------------------------------------------------------------------------
// Constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
PositionPolicy<Dimension>::
PositionPolicy():
  IncrementState<Dimension, typename Dimension::Vector>() {
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
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::position);

  // Grab the state field.
  Field<Dimension, Vector>& r = state.field(key, Vector::zero);

  // Get the FluidNodeList and check if we're enforcing compatible energy 
  // evolution or not.
  const FluidNodeList<Dimension>* nodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(r.nodeListPtr());
  CHECK(nodeListPtr != 0);

  // Get the velocity and acceleration fields.
  const KeyType velKey = State<Dimension>::buildFieldKey(HydroFieldNames::velocity, nodeListKey);
  const KeyType dvelKey = State<Dimension>::buildFieldKey(IncrementState<Dimension, Vector>::prefix() + HydroFieldNames::velocity, nodeListKey);
  CHECK(state.registered(velKey));
  CHECK(derivs.registered(dvelKey));
  const Field<Dimension, Vector>& vel = state.field(velKey, Vector::zero);
  const Field<Dimension, Vector>& dvel = derivs.field(dvelKey, Vector::zero);

  // Iterate over the internal values.
  for (int i = 0; i != nodeListPtr->numInternalNodes(); ++i) {

    // Compute time centered value for the velocity.
    const Vector vi = vel(i) + 0.5*multiplier*dvel(i);

    // Now compute the new position.
    r(i) += multiplier*vi;
   
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

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
namespace Spheral {
  template class PositionPolicy<Dim<1> >;
  template class PositionPolicy<Dim<2> >;
  template class PositionPolicy<Dim<3> >;
}
