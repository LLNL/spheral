//---------------------------------Spheral++----------------------------------//
// ReplaceState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
//
// This version assumes there is a derivative based update available.
//
// Created by JMO, Thu Aug 26 16:30:02 2004
//----------------------------------------------------------------------------//
#include "IncrementState.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/Field.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
ReplaceState(std::initializer_list<std::string> depends):
  PureReplaceState<Dimension, ValueType>(depends) {
}

//------------------------------------------------------------------------------
// Update the field as an increment
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
void
ReplaceState<Dimension, ValueType>::
updateAsIncrement(const KeyType& key,
                  State<Dimension>& state,
                  StateDerivatives<Dimension>& derivs,
                  const double multiplier,
                  const double t,
                  const double dt) {
  IncrementState<Dimension, ValueType> mIncrementStatePolicy;
  mIncrementStatePolicy.update(key, state, derivs, multiplier, t, dt);
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
bool
ReplaceState<Dimension, ValueType>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {
  return dynamic_cast<const ReplaceState<Dimension, ValueType>*>(&rhs)!= nullptr;
}

}

