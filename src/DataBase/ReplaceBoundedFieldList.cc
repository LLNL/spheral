//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//

#include "ReplaceBoundedFieldList.hh"
#include "IncrementBoundedFieldList.hh"
#include "FieldListUpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const std::string& depend0,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const std::string& depend0,
                        const std::string& depend1,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2,
                        const std::string& depend3,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2,
                        const std::string& depend3,
                        const std::string& depend4,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
ReplaceBoundedFieldList(const std::string& depend0,
                        const std::string& depend1,
                        const std::string& depend2,
                        const std::string& depend3,
                        const std::string& depend4,
                        const std::string& depend5,
                        const BoundValueType minValue,
                        const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4, depend5),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
~ReplaceBoundedFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
void
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
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
  FieldList<Dimension, ValueType> f = state.fields(fieldKey, ValueType());
  const FieldList<Dimension, ValueType> df = derivs.fields(replaceKey, ValueType());
  CHECK(f.size() == df.size());

  // Loop over the internal values of the field.
  const unsigned numNodeLists = f.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = f[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      f(k, i) = min(mMaxValue, max(mMinValue, df(k, i)));
    }
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
void
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {
  IncrementBoundedFieldList<Dimension, ValueType, BoundValueType> mIncrementFieldListPolicy(mMinValue, mMaxValue);
  mIncrementFieldListPolicy.update(key, state, derivs, multiplier, t, dt);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
bool
ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>* rhsPtr = dynamic_cast<const ReplaceBoundedFieldList<Dimension, ValueType, BoundValueType>*>(&rhs);
  return rhsPtr != 0;
}

}

