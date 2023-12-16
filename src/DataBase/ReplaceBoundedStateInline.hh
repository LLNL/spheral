//---------------------------------Spheral++----------------------------------//
// ReplaceBoundedState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//
#include "IncrementBoundedState.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"
#include "Utilities/SpheralFunctions.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(const BoundValueType minValue,
                    const BoundValueType maxValue):
  PureReplaceBoundedState<Dimension, ValueType, BoundValueType>(minValue, maxValue) {
}

template<typename Dimension, typename ValueType, typename BoundValueType>
inline
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
ReplaceBoundedState(std::initializer_list<std::string> depends,
                    const BoundValueType minValue,
                    const BoundValueType maxValue):
  PureReplaceBoundedState<Dimension, ValueType, BoundValueType>(depends, minValue, maxValue) {
}

//------------------------------------------------------------------------------
// Update the field as an increment
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType, typename BoundValueType>
inline
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
inline
bool
ReplaceBoundedState<Dimension, ValueType, BoundValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const ReplaceBoundedState<Dimension, ValueType, BoundValueType>*>(&rhs) != nullptr;
}

}
