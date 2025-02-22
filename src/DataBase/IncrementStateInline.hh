//---------------------------------Spheral++----------------------------------//
// IncrementState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Thu Aug 26 14:28:07 2004
//----------------------------------------------------------------------------//

#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(std::initializer_list<std::string> depends,
               const bool wildCardDerivs):
  FieldUpdatePolicy<Dimension, Value>(depends),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const bool wildCardDerivs):
  FieldUpdatePolicy<Dimension, Value>({}),
  mWildCardDerivs(wildCardDerivs) {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
void
IncrementState<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  // Get the state we're updating.
  auto& f = state.field(key, Value());

  // Find all the available matching derivative Field keys.
  const auto incrementKey = prefix() + fieldKey;
  const auto allkeys = derivs.keys();
  KeyType dfKey, dfNodeListKey;
  auto numDeltaFields = 0u;
  for (const auto& key: allkeys) {
    StateBase<Dimension>::splitFieldKey(key, dfKey, dfNodeListKey);
    if (dfNodeListKey == nodeListKey and
        dfKey.compare(0, incrementKey.size(), incrementKey) == 0) {
      ++numDeltaFields;

      // This delta field matches the base of increment key, so apply it.
      const auto& df = derivs.field(key, Value());
      const auto  n = f.numInternalElements();
#pragma omp parallel for
      for (auto i = 0u; i < n; ++i) {
        f(i) += multiplier*(df(i));
      }
    }
  }

  // If we're not allowing wildcard update, there should have only be one match.
  VERIFY2(mWildCardDerivs or numDeltaFields == 1,
          "IncrementState ERROR: unable to find unique match for derivative field key " << incrementKey);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
bool
IncrementState<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const IncrementState<Dimension, Value>*>(&rhs);
  return rhsPtr != nullptr;
}

//------------------------------------------------------------------------------
// Wildcard derivs attribute.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
bool
IncrementState<Dimension, Value>::
wildCardDerivs() const {
  return mWildCardDerivs;
}

template<typename Dimension, typename Value>
inline
void
IncrementState<Dimension, Value>::
wildCardDerivs(const bool val) {
  mWildCardDerivs = val;
}

}

