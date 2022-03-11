//---------------------------------Spheral++----------------------------------//
// ReplaceWithRatioPolicy -- replaces one fieldlist with the ratio of two
// fieldlists from the state.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//

#include "GSPH/Policies/ReplaceWithRatioPolicy.hh"
#include "DataBase/IncrementFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <limits.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator):
  FieldListUpdatePolicyBase<Dimension, Value>(),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator,
                       const std::string& depend0):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator,
                       const std::string& depend0,
                       const std::string& depend1):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator,
                       const std::string& depend0,
                       const std::string& depend1,
                       const std::string& depend2):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2),
  mNumerator(numerator),
  mDenomenator(denomenator){
}

template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator,
                       const std::string& depend0,
                       const std::string& depend1,
                       const std::string& depend2,
                       const std::string& depend3):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
ReplaceWithRatioPolicy(const KeyType& numerator,
                       const KeyType& denomenator,
                 const std::string& depend0,
                 const std::string& depend1,
                 const std::string& depend2,
                 const std::string& depend3,
                 const std::string& depend4,
                 const std::string& depend5):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5),
  mNumerator(numerator),
  mDenomenator(denomenator) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplaceWithRatioPolicy<Dimension, Value>::
~ReplaceWithRatioPolicy() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
ReplaceWithRatioPolicy<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double /*multiplier*/,
       const double /*t*/,
       const double /*dt*/) {
  const auto tiny = std::numeric_limits<typename Dimension::Scalar>::epsilon();
  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  CHECK(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());

  // Find the matching replacement FieldList from the StateDerivatives.
  FieldList<Dimension, Value> f = state.fields(fieldKey, Value());
  const FieldList<Dimension, Value> numer = state.fields(mNumerator, Value());
  const FieldList<Dimension, Value> denom = state.fields(mDenomenator, Value());
  CHECK(numer.size() == denom.size());

  // Loop over the internal values of the field.
  const unsigned numNodeLists = f.size();
  for (unsigned k = 0; k != numNodeLists; ++k) {
    const unsigned n = f[k]->numInternalElements();
    for (unsigned i = 0; i != n; ++i) {
      f(k, i) = numer(k, i)/std::max(denom(k, i),tiny);
    }
  }
}


//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
ReplaceWithRatioPolicy<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an replace operator.
  const ReplaceWithRatioPolicy<Dimension, Value>* rhsPtr = dynamic_cast<const ReplaceWithRatioPolicy<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}

}

