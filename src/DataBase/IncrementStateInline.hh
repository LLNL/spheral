//---------------------------------Spheral++----------------------------------//
// IncrementState -- An implementation of UpdatePolicyBase appropriate for
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
template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState():
  FieldUpdatePolicyBase<Dimension, Value>() {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const std::string& depend0):
  FieldUpdatePolicyBase<Dimension, Value>(depend0) {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const std::string& depend0,
               const std::string& depend1):
  FieldUpdatePolicyBase<Dimension, Value>(depend0, depend1) {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const std::string& depend0,
               const std::string& depend1,
               const std::string& depend2):
  FieldUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2) {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const std::string& depend0,
               const std::string& depend1,
               const std::string& depend2,
               const std::string& depend3):
  FieldUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3) {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const std::string& depend0,
               const std::string& depend1,
               const std::string& depend2,
               const std::string& depend3,
               const std::string& depend4):
  FieldUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4) {
}

template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
IncrementState(const std::string& depend0,
               const std::string& depend1,
               const std::string& depend2,
               const std::string& depend3,
               const std::string& depend4,
               const std::string& depend5):
  FieldUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
IncrementState<Dimension, Value>::
~IncrementState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
void
IncrementState<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Find the matching derivative field from the StateDerivatives.
  KeyType incrementKey = prefix() + key;
  Field<Dimension, Value>& f = state.field(key, Value());
  const Field<Dimension, Value>& df = derivs.field(incrementKey, Value());

  // Loop over the internal values of the field.
  for (unsigned i = 0; i != f.nodeList().numInternalNodes(); ++i) {
    f(i) += multiplier*(df(i));
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
inline
bool
IncrementState<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementState<Dimension, Value>* rhsPtr = dynamic_cast<const IncrementState<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

