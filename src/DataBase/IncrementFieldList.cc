//---------------------------------Spheral++----------------------------------//
// IncrementFieldList -- An implementation of UpdatePolicyBase appropriate for
// when 'ya just want to increment by derivatives:  x1 = x0 + A*dx/dt
//
// Created by JMO, Sun Oct 27 11:32:51 PDT 2013
//----------------------------------------------------------------------------//
#include "IncrementFieldList.hh"
#include "State.hh"
#include "StateDerivatives.hh"
#include "Field/FieldList.hh"
#include "Utilities/DBC.hh"

#include <regex>
#include <vector>
using std::vector;
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Constructors.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4),
  mWildCardDerivs(wildCardDerivs) {
}

template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
IncrementFieldList(const std::string& depend0,
                   const std::string& depend1,
                   const std::string& depend2,
                   const std::string& depend3,
                   const std::string& depend4,
                   const std::string& depend5,
                   const bool wildCardDerivs):
  FieldListUpdatePolicyBase<Dimension, Value>(depend0, depend1, depend2, depend3, depend4, depend5),
  mWildCardDerivs(wildCardDerivs) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
IncrementFieldList<Dimension, Value>::
~IncrementFieldList() {
}

//------------------------------------------------------------------------------
// Update the field.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
void
IncrementFieldList<Dimension, Value>::
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

  // Get the state we're updating.
  auto f = state.fields(fieldKey, Value());
  const auto numNodeLists = f.size();

  // Find all the available matching derivative FieldList keys.
  const auto incrementKey = prefix() + fieldKey;
  // cerr << "IncrementFieldList: [" << fieldKey << "] [" << incrementKey << "] : " << endl;
  const auto allkeys = derivs.fieldKeys();
  vector<string> incrementKeys;
  for (const auto& key: allkeys) {
    // if (std::regex_search(key, std::regex("^" + incrementKey))) {
    if (key.compare(0, incrementKey.size(), incrementKey) == 0) {
      incrementKeys.push_back(key);
    }
  }
  CHECK(not incrementKeys.empty());

  // If we're not allowing wildcard update, there should only be one match.
  VERIFY2(mWildCardDerivs or incrementKeys.size() == 1,
          "IncrementFieldList ERROR: unable to find unique match for derivative field key " << incrementKey);

  // Update by each of our derivative fields.
  for (const auto& key: incrementKeys) {
    const auto df = derivs.fields(key, Value());
    CHECK(df.size() == f.size());
    for (auto k = 0u; k != numNodeLists; ++k) {
      const auto n = f[k]->numInternalElements();
      for (auto i = 0u; i != n; ++i) {
        f(k, i) += multiplier*(df(k, i));
      }
    }
  }
}

//------------------------------------------------------------------------------
// Equivalence operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
IncrementFieldList<Dimension, Value>::
operator==(const UpdatePolicyBase<Dimension>& rhs) const {

  // We're only equal if the other guy is also an increment operator.
  const IncrementFieldList<Dimension, Value>* rhsPtr = dynamic_cast<const IncrementFieldList<Dimension, Value>*>(&rhs);
  return rhsPtr != 0 and mWildCardDerivs == rhsPtr->mWildCardDerivs;
}

//------------------------------------------------------------------------------
// Wildcard derivs attribute.
//------------------------------------------------------------------------------
template<typename Dimension, typename Value>
bool
IncrementFieldList<Dimension, Value>::
wildCardDerivs() const {
  return mWildCardDerivs;
}

template<typename Dimension, typename Value>
void
IncrementFieldList<Dimension, Value>::
wildCardDerivs(const bool val) {
  mWildCardDerivs = val;
}


template class IncrementFieldList<Dim<1>, Dim<3>::Vector>;
template class IncrementFieldList<Dim<2>, Dim<3>::Vector>;

}

