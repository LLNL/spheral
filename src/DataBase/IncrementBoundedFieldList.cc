//---------------------------------Spheral++----------------------------------//
// IncrementBoundedFieldList -- An implementation of FieldListUpdatePolicyBase 
// appropriate for when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt.
// This version enforces min/max bounds on the result, and therefore is only 
// appropriate for Scalar & SymTensor values.
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "IncrementBoundedFieldList.hh"
#include "UpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

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
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const std::string& depend0,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const std::string& depend0,
                          const std::string& depend1,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const std::string& depend0,
                          const std::string& depend1,
                          const std::string& depend2,
                          const BoundValueType minValue,
                          const BoundValueType maxValue):
  FieldListUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const std::string& depend0,
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
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const std::string& depend0,
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
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
IncrementBoundedFieldList(const std::string& depend0,
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
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
~IncrementBoundedFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
void
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  CHECK(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the matching derivative FieldList from the StateDerivatives.
  KeyType incrementKey = prefix() + fieldKey;
  FieldList<Dimension, ValueType> f = state.fields(fieldKey, ValueType());
  const FieldList<Dimension, ValueType> df = derivs.fields(incrementKey, ValueType());
  CHECK(f.size() == df.size());

  // Loop over the internal values of the field.
  const unsigned numNodeLists = f.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = f[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      f(k, i) = min(mMaxValue, max(mMinValue, f(k, i) + multiplier*(df(k, i))));
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
bool
IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>* rhsPtr = dynamic_cast<const IncrementBoundedFieldList<Dimension, ValueType, BoundValueType>*>(&rhs);
  if (rhsPtr == 0) return false;

  // Ok, now do we agree on min & max?
  return (minValue() == rhsPtr->minValue()) && (maxValue() == rhsPtr->maxValue());
}

}

