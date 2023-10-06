//---------------------------------Spheral++----------------------------------//
// ReplaceAndIncrementPairFieldList -- Update policy which first replaces the
//                                     pairFieldList in question then increments
//                                     it. Naturally two derivatives fields
//                                     are required one for each of the two
//                                     steps.
//
// J.M. Pearl 2022
//----------------------------------------------------------------------------//
#include "DEM/ReplaceAndIncrementPairFieldList.hh"
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
ReplaceAndIncrementPairFieldList<Dimension, Value>::
ReplaceAndIncrementPairFieldList():
  UpdatePolicyBase<Dimension>() {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementPairFieldList<Dimension, Value>::
ReplaceAndIncrementPairFieldList(const std::string& depend0):
  UpdatePolicyBase<Dimension>(depend0 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementPairFieldList<Dimension, Value>::
ReplaceAndIncrementPairFieldList(const std::string& depend0,
                                 const std::string& depend1):
  UpdatePolicyBase<Dimension>(depend0, depend1 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementPairFieldList<Dimension, Value>::
ReplaceAndIncrementPairFieldList(const std::string& depend0,
                                 const std::string& depend1,
                                 const std::string& depend2):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementPairFieldList<Dimension, Value>::
ReplaceAndIncrementPairFieldList(const std::string& depend0,
                                 const std::string& depend1,
                                 const std::string& depend2,
                                 const std::string& depend3):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementPairFieldList<Dimension, Value>::
ReplaceAndIncrementPairFieldList(const std::string& depend0,
                                 const std::string& depend1,
                                 const std::string& depend2,
                                 const std::string& depend3,
                                 const std::string& depend4):
  UpdatePolicyBase<Dimension>(depend0, depend1, depend2, depend3, depend4 ) {
}

template<typename Dimension, typename Value>
ReplaceAndIncrementPairFieldList<Dimension, Value>::
ReplaceAndIncrementPairFieldList(const std::string& depend0,
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
ReplaceAndIncrementPairFieldList<Dimension, Value>::
~ReplaceAndIncrementPairFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
ReplaceAndIncrementPairFieldList<Dimension, Value>::
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
ReplaceAndIncrementPairFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const ReplaceAndIncrementPairFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const ReplaceAndIncrementPairFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0;
}


}

