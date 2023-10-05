//---------------------------------Spheral++----------------------------------//
// ReplaceState -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to replace the state value with the new.
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
ReplaceState():
  FieldUpdatePolicy<Dimension>() {
}

template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0):
  FieldUpdatePolicy<Dimension>(depend0) {
}

template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1):
  FieldUpdatePolicy<Dimension>(depend0, depend1) {
}

template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2) {
}

template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2,
             const std::string& depend3):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3) {
}

template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2,
             const std::string& depend3,
             const std::string& depend4):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3, depend4) {
}

template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
ReplaceState(const std::string& depend0,
             const std::string& depend1,
             const std::string& depend2,
             const std::string& depend3,
             const std::string& depend4,
             const std::string& depend5):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
ReplaceState<Dimension, ValueType>::
~ReplaceState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename ValueType>
inline
void
ReplaceState<Dimension, ValueType>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  // Find the matching replacement field from the StateDerivatives.
  const auto  replaceKey = prefix() + key;
  auto&       f = state.field(key, ValueType());
  const auto& df = derivs.field(replaceKey, ValueType());

  // Loop over the internal values of the field.
  const auto n = f.nodeList().numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    f(i) = df(i);
  }
}

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

  // We're only equal if the other guy is also an increment operator.
  const auto rhsPtr = dynamic_cast<const ReplaceState<Dimension, ValueType>*>(&rhs);
  return rhsPtr != nullptr;
}

}

