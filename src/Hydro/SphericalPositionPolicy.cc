//---------------------------------Spheral++----------------------------------//
// SphericalPositionPolicy
//
// Specializes the IncrementFieldListPolicy for use updating the position in
// spherical coordinates.  This position does not allow points to pass through
// the origin.
//
// Created by JMO, Tue Mar 29 11:12:31 PDT 2022
//----------------------------------------------------------------------------//
#include "SphericalPositionPolicy.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <regex>
#include <vector>
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructor
//------------------------------------------------------------------------------
SphericalPositionPolicy::
SphericalPositionPolicy():
  UpdatePolicyBase<Dim<1>>() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
SphericalPositionPolicy::
~SphericalPositionPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
void
SphericalPositionPolicy::
update(const KeyType& key,
       State<Dim<1>>& state,
       StateDerivatives<Dim<1>>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Get the state we're updating.
  auto f = state.fields(fieldKey, Vector::zero);
  const auto numNodeLists = f.size();

  // Find all the available matching derivative Field keys.
  const auto incrementKey = prefix() + fieldKey;
  const auto allkeys = derivs.fullFieldKeys();
  vector<string> incrementKeys;
  for (const auto& key: allkeys) {
    if (key.compare(0, incrementKey.size(), incrementKey) == 0) {
      incrementKeys.push_back(key);
    }
  }
  CHECK(not incrementKeys.empty());

  // Update by each of our derivative fields.
  for (const auto& key: incrementKeys) {
    const auto df = derivs.fields(key, Vector::zero);
    CHECK(df.size() == f.size());
    for (auto k = 0u; k != numNodeLists; ++k) {
      const auto n = f[k]->numInternalElements();
      for (auto i = 0u; i != n; ++i) {
        // This is where we diverge from the standard IncrementState.  Ensure we cannot cross to
        // negative radius.
        f(k,i) = std::max(0.5*f(k,i), f(k,i) + multiplier*(df(k, i)));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
bool
SphericalPositionPolicy::
operator==(const UpdatePolicyBase<Dim<1>>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const SphericalPositionPolicy*>(&rhs);
  return (rhsPtr != nullptr);
}

}

