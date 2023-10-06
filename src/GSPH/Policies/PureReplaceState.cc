//---------------------------------Spheral++----------------------------------//
// PureReplaceState -- replaces one fieldlists values with those of 
//                         another specified by its key
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/Policies/PureReplaceState.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey):
  ReplaceState<Dimension, Value>(),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0):
  ReplaceState<Dimension, Value>(depend0),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1):
  ReplaceState<Dimension, Value>(depend0, depend1),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2):
  ReplaceState<Dimension, Value>(depend0, depend1, depend2),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceState<Dimension, Value>::
PureReplaceState(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3):
  ReplaceState<Dimension, Value>(depend0, depend1, depend2, depend3),
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
  ReplaceState<Dimension, Value>(depend0, depend1, depend2, depend3, depend4),
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
  ReplaceState<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5),
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

