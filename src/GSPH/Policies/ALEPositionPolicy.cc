//---------------------------------Spheral++----------------------------------//
// ALEPositionPolicy -- This is basically a direct copy of the standard 
//                      position policy but instead we're substituting in 
//                      the nodal velocity as the derivative.
//
// J. M. Pearl 2023
//----------------------------------------------------------------------------//

#include "ALEPositionPolicy.hh"
#include "GSPH/GSPHFieldNames.hh"
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
ALEPositionPolicy<Dimension>::
ALEPositionPolicy():
  IncrementFieldList<Dimension, typename Dimension::Vector>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
ALEPositionPolicy<Dimension>::
~ALEPositionPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
ALEPositionPolicy<Dimension>::
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
  const auto vel = state.fields(GSPHFieldNames::nodalVelocity, Vector::zero);

  std::cout << "ALE POSITION UPDATE" << std::endl;
  // Walk the fields.
  for (auto i = 0u; i != numFields; ++i) {
    const auto n = r[i]->numInternalElements();
    for (auto j = 0u; j < n; ++j) {
      r(i,j) += multiplier*vel(i,j);
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
ALEPositionPolicy<Dimension>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ALEPositionPolicy<Dimension>* rhsPtr = dynamic_cast<const ALEPositionPolicy<Dimension>*>(&rhs);
  if (rhsPtr == 0) {
    return false;
  } else {
    return true;
  }
}

}

