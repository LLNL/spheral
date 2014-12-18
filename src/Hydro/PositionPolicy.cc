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

using namespace std;
using FieldSpace::FieldList;
using NodeSpace::NodeList;
using NodeSpace::FluidNodeList;

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
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(fieldKey == HydroFieldNames::position and 
          nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  FieldList<Dimension, Vector> r = state.fields(fieldKey, Vector::zero);
  const unsigned numFields = r.numFields();

  // Get the velocity and acceleration fields.
  const FieldList<Dimension, Vector> vel = state.fields(HydroFieldNames::velocity, Vector::zero);
  const FieldList<Dimension, Vector> dvel = derivs.fields(IncrementFieldList<Dimension, Vector>::prefix() + HydroFieldNames::velocity, Vector::zero);

  // Walk the fields.
  for (unsigned i = 0; i != numFields; ++i) {

    // Get the FluidNodeList and check if we're enforcing compatible energy 
    // evolution or not.
    const FluidNodeList<Dimension>* nodeListPtr = dynamic_cast<const FluidNodeList<Dimension>*>(r[i]->nodeListPtr());
    CHECK(nodeListPtr != 0);

    // Iterate over the internal values.
    for (unsigned j = 0; j != r[i]->numInternalElements(); ++j) {

      // Compute time centered value for the velocity.
      const Vector vi = vel(i,j) + 0.5*multiplier*dvel(i,j);

      // Now compute the new position.
      r(i,j) += multiplier*vi;
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

