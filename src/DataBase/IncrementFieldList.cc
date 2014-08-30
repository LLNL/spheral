//---------------------------------Spheral++----------------------------------//
// IncrementFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "IncrementFieldList.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

using FieldSpace::FieldList;

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList():
  FieldListUpdatePolicyBase<Dimension, Value>() {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
               const std::string& depend1):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
               const std::string& depend1,
               const std::string& depend2):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
               const std::string& depend1,
               const std::string& depend2,
               const std::string& depend3):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
               const std::string& depend1,
               const std::string& depend2,
               const std::string& depend3,
               const std::string& depend4):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
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
IncrementFieldList<Dimension, Value>::
~IncrementFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
IncrementFieldList<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the matching derivative FieldList from the StateDerivatives.
  KeyType incrementKey = prefix() + fieldKey;
  FieldSpace::FieldList<Dimension, Value> f = state.fields(fieldKey, Value());
  const FieldSpace::FieldList<Dimension, Value> df = derivs.fields(incrementKey, Value());
  CHECK(f.size() == df.size());

  // Loop over the internal values of the field.
  const unsigned numNodeLists = f.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = f[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      f(k, i) += multiplier*(df(k, i));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
IncrementFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const IncrementFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

