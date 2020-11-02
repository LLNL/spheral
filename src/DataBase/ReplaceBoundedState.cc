//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//
#include "ReplaceBoundedState.hh"
#include "IncrementBoundedState.hh"
#include "FieldUpdatePolicyBase.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
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
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const std::string& depend3,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const std::string& depend3,
                    const std::string& depend4,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const std::string& depend0,
                    const std::string& depend1,
                    const std::string& depend2,
                    const std::string& depend3,
                    const std::string& depend4,
                    const std::string& depend5,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  FieldUpdatePolicyBase<Dimension, ValueType>(depend0, depend1, depend2, depend3, depend4, depend5),
  mMinValue(minValue),
  mMaxValue(maxValue) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
~ReplaceBoundedState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
void
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  // Find the matching replacement field from the StateDerivatives.
  KeyType replaceKey = prefix() + key;
  Field<Dimension, ValueType>& f = state.field(key, ValueType());
  const Field<Dimension, ValueType>& df = derivs.field(replaceKey, ValueType());

  // Loop over the internal values of the field.
  for (auto i = 0u; i != f.nodeList().numInternalNodes(); ++i) {
    f(i) = min(mMaxValue, max(mMinValue, df(i)));
  }
}

template<typename Dimension, typename ValueType, typename BoundValueType>
void
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {
  IncrementBoundedState<Dimension, ValueType, BoundValueType> mIncrementStatePolicy(mMinValue, mMaxValue);
  mIncrementStatePolicy.update(key, state, derivs, multiplier, t, dt);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
bool
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ReplaceBoundedState<Dimension, ValueType, BoundValueType>* rhsPtr = dynamic_cast<const ReplaceBoundedState<Dimension, ValueType, BoundValueType>*>(&rhs);
  return rhsPtr != 0;
}

}

