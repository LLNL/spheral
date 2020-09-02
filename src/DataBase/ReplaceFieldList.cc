//---------------------------------Spheral++----------------------------------//
// ReplaceFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "ReplaceFieldList.hh"
#include "IncrementFieldList.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList():
  FieldListUpdatePolicyBase<Dimension, Value>() {
}

template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const std::string& depend0):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0) {
}

template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const std::string& depend0,
                 const std::string& depend1):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1) {
}

template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2) {
}

template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3) {
}

template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4) {
}

template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
ReplaceFieldList(const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4,
                 const std::string& depend5):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplaceFieldList<Dimension, Value>::
~ReplaceFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
ReplaceFieldList<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  CHECK(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the matching replacement FieldList from the StateDerivatives.
  KeyType replaceKey = prefix() + fieldKey;
  FieldList<Dimension, Value> f = state.fields(fieldKey, Value());
  const FieldList<Dimension, Value> df = derivs.fields(replaceKey, Value());
  CHECK(f.size() == df.size());

  // Loop over the internal values of the field.
  const unsigned numNodeLists = f.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = f[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      f(k, i) = df(k, i);
    }
  }
}

template<typename Dimension, typename Value>
void
ReplaceFieldList<Dimension, Value>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {
  IncrementFieldList<Dimension, Value> mIncrementStatePolicy;
  mIncrementStatePolicy.update(key, state, derivs, multiplier, t, dt);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
ReplaceFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const ReplaceFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const ReplaceFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

