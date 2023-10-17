//---------------------------------Spheral++----------------------------------//
// PureReplaceState -- replaces one fieldlists values with those of 
//                         another specified by its key
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/Policies/PureReplaceState.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey):
  FieldUpdatePolicy<Dimension>(),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0):
  FieldUpdatePolicy<Dimension>(depend0),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1):
  FieldUpdatePolicy<Dimension>(depend0, depend1),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3, depend4),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4,
                 const std::string& depend5):
  FieldUpdatePolicy<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5),
  mReplaceKey(derivFieldListKey) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
~PureReplaceState() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
PureReplaceState<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {

  // Split the key into Field and NodeList keys
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);

  // The state we're updating
  auto& f = state.field(key, Value());

  // Find the matching replacement field from the StateDerivatives.
  const auto  buildKey = [&](const std::string& fkey) { return StateBase<Dimension>::buildFieldKey(fkey, nodeListKey); };
  const auto& df = derivs.field(buildKey(mReplaceKey), Value());

  // Loop over the internal values of the field.
  const auto n = f.nodeList().numInternalNodes();
#pragma omp parallel for
  for (auto i = 0u; i < n; ++i) {
    f(i) = df(i);
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
PureReplaceState<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const auto* rhsPtr = dynamic_cast<const PureReplaceState<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

