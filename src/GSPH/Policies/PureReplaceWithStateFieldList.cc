//---------------------------------Spheral++----------------------------------//
// PureReplaceWithStateFieldList -- replaces one fieldlists values with  
//                                  a state field list specified by its key
//
// J.M. Pearl 2023
//----------------------------------------------------------------------------//

#include "GSPH/Policies/PureReplaceWithStateFieldList.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
PureReplaceWithStateFieldList(const KeyType& derivFieldListKey):
  FieldListUpdatePolicyBase<Dimension, Value>(),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
PureReplaceWithStateFieldList(const KeyType& derivFieldListKey,
                     const std::string& depend0):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
PureReplaceWithStateFieldList(const KeyType& derivFieldListKey,
                     const std::string& depend0,
                     const std::string& depend1):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
PureReplaceWithStateFieldList(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
PureReplaceWithStateFieldList(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
PureReplaceWithStateFieldList(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4),
  mReplaceKey(derivFieldListKey) {
}

template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
PureReplaceWithStateFieldList(const KeyType& derivFieldListKey,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4,
                 const std::string& depend5):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5),
  mReplaceKey(derivFieldListKey) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
PureReplaceWithStateFieldList<Dimension, Value>::
~PureReplaceWithStateFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
PureReplaceWithStateFieldList<Dimension, Value>::
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
  FieldList<Dimension, Value> f = state.fields(fieldKey, Value());

  // field we're replacing it with
  const FieldList<Dimension, Value> df = state.fields(mReplaceKey, Value());
  CHECK(f.size() == df.size());

  // Loop over the internal values of the field.
  const unsigned numNodeLists = f.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = f[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      f(k, i) = df(k, i);
    }
  }
}


//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
PureReplaceWithStateFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const PureReplaceWithStateFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const PureReplaceWithStateFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

