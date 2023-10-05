//---------------------------------Spheral++----------------------------------//
// IncrementBoundedState -- An implementation of UpdatePolicyBase appropriate for
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
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension>(),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension>(depend0),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const std::string& depend1,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension>(depend0, depend1),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const std::string& depend1,
                      const std::string& depend2,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const std::string& depend1,
                      const std::string& depend2,
                      const std::string& depend3,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const std::string& depend1,
                      const std::string& depend2,
                      const std::string& depend3,
                      const std::string& depend4,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3, depend4),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
IncrementBoundedState(const std::string& depend0,
                      const std::string& depend1,
                      const std::string& depend2,
                      const std::string& depend3,
                      const std::string& depend4,
                      const std::string& depend5,
                      const BoundValueType minValue,
                      const BoundValueType maxValue):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
~IncrementBoundedState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
void
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Find the matching derivative field from the StateDerivatives.
  const auto  incrementKey = prefix() + key;
  auto&       f = state.field(key, ValueType());
  const auto& df = derivs.field(incrementKey, ValueType());

  // Loop over the internal values of the field.
  const auto n = f.nodeList().numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    f(i) = std::min(mMaxValue, std::max(mMinValue, f(i) + multiplier*(df(i))));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
bool
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const IncrementBoundedState<Dimension, ValueType, BoundValueType>*>(&rhs);
  if (rhsPtr == nullptr) return false;

  // Ok, now do we agree on min & max?
  return (minValue() == rhsPtr->minValue()) && (maxValue() == rhsPtr->maxValue());
}

//------------------------------------------------------------------------------
// Min value.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
BoundValueType
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
minValue() const {
  return mMinValue;
}

//------------------------------------------------------------------------------
// Max value.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
BoundValueType
IncrementBoundedState<Dimension, ValueType, BoundValueType>::
maxValue() const {
  return mMaxValue;
}

}
