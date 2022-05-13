//---------------------------------Spheral++----------------------------------//
// ReplaceAndIncrementFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//----------------------------------------------------------------------------//
#include "DEM/ReplaceAndIncrementFieldList.hh"
#include "DataBase/State.hh"
#include "DataBase/StateDerivatives.hh"
#include "Field/FieldList.hh"

#include <vector>
using std::vector;
using std::string;


namespace Spheral {

//------------------------------------------------------------------------------
// Constructors mostly wrappers to the base
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
ReplaceAndIncrementFieldList():
  FieldListUpdatePolicyBase<Dimension, Value>() {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
ReplaceAndIncrementFieldList(const std::string& depend0):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
ReplaceAndIncrementFieldList(const std::string& depend0,
                   const std::string& depend1):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
ReplaceAndIncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
ReplaceAndIncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
ReplaceAndIncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
ReplaceAndIncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4,
                   const std::string& depend5):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5 ) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
ReplaceAndIncrementFieldList<Dimension, Value>::
~ReplaceAndIncrementFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
ReplaceAndIncrementFieldList<Dimension, Value>::
update(const KeyType& key,
       State<Dimension>& state,
       StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double /*t*/,
       const double /*dt*/) {

  // Get the field name portion of the key.
  KeyType fieldKey, nodeListKey;
  StateBase<Dimension>::splitFieldKey(key, fieldKey, nodeListKey);
  REQUIRE(nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
  
  // get our replace field
  KeyType incrementKey = incrementPrefix() + fieldKey;
  KeyType replaceKey = replacePrefix() + fieldKey;
  auto f = state.fields(fieldKey, Value());
  const auto newf = derivs.fields(replaceKey, Value());
  const auto df = derivs.fields(incrementKey, Value());
  
  CHECK(f.size() == df.size());
  CHECK(f.size() == newf.size());

  const auto numNodeLists = f.size();

  // Update 
  for (auto k = 0u; k != numNodeLists; ++k){
    const auto numNodes = f[k]->numInternalElements();
    for (auto i = 0u; i != numNodes; ++i){
      const auto numContacts = df(k,i).size();
      for (auto j = 0u; j != numContacts; ++j){
        f(k, i)[j] = newf(k, i)[j] + multiplier*(df(k, i)[j]);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
ReplaceAndIncrementFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ReplaceAndIncrementFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const ReplaceAndIncrementFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}


}

