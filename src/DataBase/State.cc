//---------------------------------Spheral++----------------------------------//
// State -- Accumulate and cart the state for a set of physics packages around.
//
// Created by JMO, Fri Aug 27 10:56:40 2004
//----------------------------------------------------------------------------//
#include "State.hh"
#include "StateBase.hh"
#include "StateDerivatives.hh"
#include "DataBase.hh"
#include "UpdatePolicyBase.hh"
#include "FieldUpdatePolicyBase.hh"
#include "Physics/Physics.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Geometry/Dimension.hh"

#include <string>
using std::vector;
using std::map;
using std::set;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

namespace {
//------------------------------------------------------------------------------
// Helper to compare a set of state to a key.
// We assume that wildcards are only applied at the end of a Key!
//------------------------------------------------------------------------------
// template<typename Key, typename Iterator>
// bool keyPresent(const Key& key,
//                 const Key& wildcard,
//                 const Iterator setBegin,
//                 const Iterator setEnd) {
//   const size_t pos = key.find(wildcard);
//   if (pos == string::npos) {
//     Key fieldKey, nodeListKey;
//     Iterator itr = setBegin;
//     while (itr != setEnd) {
//       StateBase<Dim<1> >::splitFieldKey(*itr, fieldKey, nodeListKey);
//       if (fieldKey == key) return true;
//       ++itr;
//     }
//     return false;
//     //return (find(setBegin, setEnd, key) != setEnd);
//   } else {
//     const Key testKey = key.substr(0, pos);
//     Iterator itr = setBegin;
//     while (itr != setEnd) {
//       if (testKey == itr->substr(0, pos)) return true;
//       ++itr;
//     }
//     return false;
//   }
// }

// //------------------------------------------------------------------------------
// // Helper to compare a set of state to a key.
// // We assume that wildcards are only applied at the end of a Key!
// //------------------------------------------------------------------------------
// bool compareKeysByField(const std::string& a,
//                         const std::string& b) {
//   std::string af, bf, an, bn;
//   StateBase<Dim<1> >::splitFieldKey(a, af, an);
//   StateBase<Dim<1> >::splitFieldKey(b, bf, bn);
//   return af < bf;
// }

}

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>::
State():
  StateBase<Dimension>(),
  mPolicyMap(),
  mTimeAdvanceOnly(false) {
}

//------------------------------------------------------------------------------
// Construct with the state for the given set of Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>::
State(DataBase<Dimension>& dataBase,
      typename State<Dimension>::PackageList& physicsPackages):
  StateBase<Dimension>(),
  mPolicyMap(),
  mTimeAdvanceOnly(false) {
  // Iterate over the physics packages, and have them register their state.
  for (PackageIterator itr = physicsPackages.begin();
       itr != physicsPackages.end();
       ++itr) (*itr)->registerState(dataBase, *this);
}

//------------------------------------------------------------------------------
// Construct with the state for the given set of Physics packages.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>::
State(DataBase<Dimension>& dataBase,
      typename State<Dimension>::PackageIterator physicsPackageBegin,
      typename State<Dimension>::PackageIterator physicsPackageEnd):
  StateBase<Dimension>(),
  mPolicyMap(),
  mTimeAdvanceOnly(false) {
  // Iterate over the physics packages, and have them register their state.
  for (PackageIterator itr = physicsPackageBegin;
       itr != physicsPackageEnd;
       ++itr) (*itr)->registerState(dataBase, *this);
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>::
State(const State<Dimension>& rhs):
  StateBase<Dimension>(rhs),
  mPolicyMap(rhs.mPolicyMap),
  mTimeAdvanceOnly(rhs.mTimeAdvanceOnly) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>::
~State() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
State<Dimension>&
State<Dimension>::
operator=(const State<Dimension>& rhs) {
  if (this != &rhs) {
    StateBase<Dimension>::operator=(rhs);
    mPolicyMap = rhs.mPolicyMap;
    mTimeAdvanceOnly = rhs.mTimeAdvanceOnly;
  }
  return *this;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension>
bool
State<Dimension>::
operator==(const StateBase<Dimension>& rhs) const {
  return StateBase<Dimension>::operator==(rhs);
}

//------------------------------------------------------------------------------
// Update the state with the given derivatives object, according to the per
// state field policies.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
update(StateDerivatives<Dimension>& derivs,
       const double multiplier,
       const double t,
       const double dt) {

  // cerr << "################################################################################" << endl;

  // Prepare lists of the keys to be completed.
  vector<KeyType> fieldsToBeCompleted;
  map<KeyType, set<KeyType> > stateToBeCompleted;
  for (typename PolicyMapType::const_iterator itr = mPolicyMap.begin();
       itr != mPolicyMap.end();
       ++itr) {
    for (typename map<KeyType, PolicyPointer>::const_iterator pitr = itr->second.begin();
         pitr != itr->second.end();
         ++pitr) {
      stateToBeCompleted[itr->first].insert(pitr->first);
    }
  }
  for (typename map<KeyType, set<KeyType> >::const_iterator itr = stateToBeCompleted.begin();
       itr != stateToBeCompleted.end();
       ++itr) fieldsToBeCompleted.push_back(itr->first);
  CHECK(fieldsToBeCompleted.size() == stateToBeCompleted.size());

  // Iterate until all state has been updated.
  while (not stateToBeCompleted.empty()) {

    // Walk the remaining state to be completed.
    vector<KeyType> stateToRemove;
    for (typename map<KeyType, set<KeyType> >::iterator itr = stateToBeCompleted.begin();
         itr != stateToBeCompleted.end();
         ++itr) {
      const KeyType fieldKey = itr->first;
      const set<KeyType>& remainingKeys = stateToBeCompleted[fieldKey];

      // Walk the remaining individual keys for this fieldKey.
      for (typename set<KeyType>::const_iterator itr = remainingKeys.begin();
           itr != remainingKeys.end();
           ++itr) {
        const KeyType key = *itr;
        const PolicyPointer policyPtr = mPolicyMap[fieldKey][key];

        // Check if all the dependencies for this state have been satisfied yet.
        vector<KeyType> unmetDependencies;
        set_intersection(fieldsToBeCompleted.begin(), fieldsToBeCompleted.end(),
                         policyPtr->dependencies().begin(), policyPtr->dependencies().end(),
                         back_inserter(unmetDependencies));
        if (unmetDependencies.empty()) {

          // We also require that any FieldList policies fire before NodeList specific 
          // versions of the same Field names.
          KeyType fieldKey, nodeListKey;
          this->splitFieldKey(key, fieldKey, nodeListKey);
          bool fire = (nodeListKey == UpdatePolicyBase<Dimension>::wildcard());
          if (not fire) {
            const KeyType wildKey = this->buildFieldKey(fieldKey, UpdatePolicyBase<Dimension>::wildcard());
            fire = (find(remainingKeys.begin(), remainingKeys.end(), wildKey) == remainingKeys.end());
          }

          if (fire) {
            // cerr <<" --> Update " << key << endl;
            if (mTimeAdvanceOnly) {
              policyPtr->updateAsIncrement(key, *this, derivs, multiplier, t, dt);
            } else {
              policyPtr->update(key, *this, derivs, multiplier, t, dt);
            }

            // List this field as completed.
            stateToRemove.push_back(key);
          }
        }
      }
    }

    // Check that the some state has been updated on this iteration.  If not, then *somebody*
    // has specified a circular dependency tree!
    if (stateToRemove.empty()) {
      std::stringstream message;
      message << "State::update ERROR: someone has specified a circular state dependency.\n"
              << "Remaining State:\n";
      for (typename map<KeyType, set<KeyType> >::const_iterator itr = stateToBeCompleted.begin();
           itr != stateToBeCompleted.end();
           ++itr) message << "   " << itr->first << "\n";
      message << "State dependencies:\n";
      for (typename PolicyMapType::iterator itr = mPolicyMap.begin();
           itr != mPolicyMap.end();
           ++itr) {
        const KeyType fieldKey = itr->first;
        const map<KeyType, PolicyPointer>& keysAndPolicies = itr->second;
        for (typename map<KeyType, PolicyPointer>::const_iterator pitr = keysAndPolicies.begin();
             pitr != keysAndPolicies.end();
             ++pitr) {
          const KeyType key = pitr->first;
          const PolicyPointer policyPtr = pitr->second;
          message << key << " : ";
          for (typename vector<KeyType>::const_iterator depItr = policyPtr->dependencies().begin();
               depItr != policyPtr->dependencies().end();
               ++depItr) message << *depItr << " ";
          message << "\n";
        }
      }
      VERIFY2(not stateToRemove.empty(), message.str());
    }
    
    // Remove the completed state.
    for (typename vector<KeyType>::const_iterator keyItr = stateToRemove.begin();
         keyItr != stateToRemove.end();
         ++keyItr) {
      KeyType fieldKey, nodeListKey;
      this->splitFieldKey(*keyItr, fieldKey, nodeListKey);
      CHECK(stateToBeCompleted.find(fieldKey) != stateToBeCompleted.end());
      stateToBeCompleted[fieldKey].erase(*keyItr);
      if (stateToBeCompleted[fieldKey].empty()) stateToBeCompleted.erase(fieldKey);
    }
    fieldsToBeCompleted = vector<KeyType>();
    for (typename map<KeyType, set<KeyType> >::const_iterator itr = stateToBeCompleted.begin();
         itr != stateToBeCompleted.end();
         ++itr) fieldsToBeCompleted.push_back(itr->first);
    CHECK(fieldsToBeCompleted.size() == stateToBeCompleted.size());
  }
}

//------------------------------------------------------------------------------
// The set of keys for all registered policies.
//------------------------------------------------------------------------------
template<typename Dimension>
vector<typename State<Dimension>::KeyType>
State<Dimension>::
policyKeys() const {
  vector<KeyType> result;
  for (typename PolicyMapType::const_iterator itr = mPolicyMap.begin();
       itr != mPolicyMap.end();
       ++itr) result.push_back(itr->first);
  ENSURE(result.size() == mPolicyMap.size());
  return result;
}

//------------------------------------------------------------------------------
// Return the policy for the given key.
//------------------------------------------------------------------------------
template<typename Dimension>
typename State<Dimension>::PolicyPointer
State<Dimension>::
policy(const typename State<Dimension>::KeyType& key) const {
  KeyType fieldKey, nodeKey;
  this->splitFieldKey(key, fieldKey, nodeKey);
  const typename PolicyMapType::const_iterator outerItr = mPolicyMap.find(fieldKey);
  if (outerItr == mPolicyMap.end()) return PolicyPointer();
  // VERIFY2(outerItr != mPolicyMap.end(),
  //         "State ERROR: attempted to retrieve non-existent policy for key " << key);
  const std::map<KeyType, PolicyPointer>& policies = outerItr->second;
  const typename std::map<KeyType, PolicyPointer>::const_iterator innerItr = policies.find(key);
  if (innerItr == policies.end()) return PolicyPointer();
  // VERIFY2(innerItr != policies.end(),
  //         "State ERROR: attempted to retrieve non-existent policy for key " << key);
  return innerItr->second;
}

//------------------------------------------------------------------------------
// Remove the policy associated with the given key.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
removePolicy(const typename State<Dimension>::KeyType& key) {
  KeyType fieldKey, nodeKey;
  this->splitFieldKey(key, fieldKey, nodeKey);
  typename PolicyMapType::iterator outerItr = mPolicyMap.find(fieldKey);
  VERIFY2(outerItr != mPolicyMap.end(),
          "State ERROR: attempted to remove non-existent policy for field key " << fieldKey);
  std::map<KeyType, PolicyPointer>& policies = outerItr->second;
  typename std::map<KeyType, PolicyPointer>::iterator innerItr = policies.find(key);
  if (innerItr == policies.end()) {
    cerr << "State ERROR: attempted to remove non-existent policy for inner key " << key << endl
         << "Known keys are: " << endl;
    for (auto itr = policies.begin(); itr != policies.end(); ++itr) cerr << " --> " << itr->first << endl;
    VERIFY(innerItr != policies.end());
  }
  policies.erase(innerItr);
  if (policies.size() == 0) mPolicyMap.erase(outerItr);
}

//------------------------------------------------------------------------------
// Remove the policy associated with a Field.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
removePolicy(FieldBase<Dimension>& field) {
  this->removePolicy(StateBase<Dimension>::key(field));
}

//------------------------------------------------------------------------------
// Remove the policy associated with a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
State<Dimension>::
removePolicy(FieldListBase<Dimension>& fieldList) {
  this->removePolicy(StateBase<Dimension>::key(fieldList));
}

}

