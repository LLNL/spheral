//---------------------------------Spheral++----------------------------------//
// ReplaceAndIncrementPairFieldList -- Update policy which replaces the values
//                                     of a pairFieldList.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#include "DEM/ReplacePairFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"

using std::string;


namespace Spheral {

//------------------------------------------------------------------------------
// Constructors mostly wrappers to the base
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
ReplacePairFieldList():
  UpdatePolicyBase<Dimension>() {
}

template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
ReplacePairFieldList(const std::string& depend0):
  UpdatePolicyBase<Dimension>(depend0 ) {
}

template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
ReplacePairFieldList(const std::string& depend0,
                   const std::string& depend1):
  UpdatePolicyBase<Dimension>(depend0, depend1 ) {
}

template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
ReplacePairFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2 ) {
}

template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
ReplacePairFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3 ) {
}

template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
ReplacePairFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4 ) {
}

template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
ReplacePairFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4,
                   const std::string& depend5):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4, depend5 ) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplacePairFieldList<Dimension, Value>::
~ReplacePairFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
ReplacePairFieldList<Dimension, Value>::
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
  KeyType replaceKey = prefix() + fieldKey;
  auto f = state.fields(fieldKey, Value());
  const auto df = derivs.fields(replaceKey, Value());
  CHECK(f.size() == df.size());

  // Loop over the internal values of the field.
  const auto numNodeLists = f.size();
  for (auto k = 0u; k != numNodeLists; ++k) {
    const auto numNodes = f[k]->numInternalElements();
    for (auto i = 0u; i != numNodes; ++i) {
      const auto numContacts = df(k,i).size();
      if (f(k,i).size()!=numContacts) f(k,i).resize(numContacts);
      for (auto j = 0u; j != numContacts; ++j) {
        f(k, i)[j] = df(k, i)[j];
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
ReplacePairFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an Replace operator.
  const auto* rhsPtr = dynamic_cast<const ReplacePairFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != nullptr;
}


}

