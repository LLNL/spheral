//---------------------------------Spheral++----------------------------------//
// CompositeFieldListPolicy -- An implementation of UpdatePolicyBase which
// consists of a collection of individual Field policies that should match
// the Fields in a FieldList.
//
// Created by JMO, Sun Nov  3 14:11:32 PST 2013
//----------------------------------------------------------------------------//
#include "FieldListUpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
CompositeFieldListPolicy<Dimension, ValueType>::
CompositeFieldListPolicy():
  FieldListUpdatePolicyBase<Dimension, ValueType>(),
  mPolicyPtrs() {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
CompositeFieldListPolicy<Dimension, ValueType>::
~CompositeFieldListPolicy() {
}

//------------------------------------------------------------------------------
// Update the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
void
CompositeFieldListPolicy<Dimension, ValueType>::
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

  // Find the FieldList for this key.
  FieldList<Dimension, ValueType> f = state.fields(fieldKey, ValueType());

  // Walk the paired Fields and policies.
  const unsigned numFields = f.numFields();
  CHECK(mPolicyPtrs.size() == numFields);
  for (unsigned i = 0; i != numFields; ++i) {
    KeyType fkey = StateBase<Dimension>::buildFieldKey(fieldKey, f[i]->nodeListPtr()->name());
    mPolicyPtrs[i]->update(fkey, state, derivs, multiplier, t, dt);
  }
}

//------------------------------------------------------------------------------
// Update the FieldList as increment.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
void
CompositeFieldListPolicy<Dimension, ValueType>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the FieldList for this key.
  FieldList<Dimension, ValueType> f = state.fields(fieldKey, ValueType());

  // Walk the paired Fields and policies.
  const unsigned numFields = f.numFields();
  CHECK(mPolicyPtrs.size() == numFields);
  for (unsigned i = 0; i != numFields; ++i) {
    KeyType fkey = StateBase<Dimension>::buildFieldKey(fieldKey, f[i]->nodeListPtr()->name());
    mPolicyPtrs[i]->updateAsIncrement(fkey, state, derivs, multiplier, t, dt);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
bool
CompositeFieldListPolicy<Dimension, ValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const CompositeFieldListPolicy<Dimension, ValueType>* rhsPtr = dynamic_cast<const CompositeFieldListPolicy<Dimension, ValueType>*>(&rhs);
  if (rhsPtr == 0) return false;
  return (this->mPolicyPtrs == rhsPtr->mPolicyPtrs);
}

//------------------------------------------------------------------------------
// Add a new UpdatePolicy to this thing.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
void
CompositeFieldListPolicy<Dimension, ValueType>::
push_back(UpdatePolicyBase<Dimension>* policyPtr) {
  mPolicyPtrs.push_back(std::unique_ptr<UpdatePolicyBase<Dimension> >(policyPtr));
  const auto newdeps = policyPtr->dependencies();
  if (mPolicyPtrs.size() == 1) {
    for (auto itr = newdeps.begin();
         itr != newdeps.end();
         ++itr) this->addDependency(*itr);
  }
  ENSURE(newdeps == this->dependencies());
}

}

