//---------------------------------Spheral++----------------------------------//
// StateBase -- The base class for State and StateDerivatives method, providing
// the common methods for registering and storing Fields, as well as serving
// up FieldLists based on these Fields.
//
// Created by JMO, Wed Aug 25 22:23:35 2004
//----------------------------------------------------------------------------//
#include "StateBase.hh"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/DBC.hh"

#include <algorithm>
#include <sstream>
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::min;
using std::max;
using std::abs;

namespace Spheral {

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateBase<Dimension>::
StateBase():
  mStorage(),
  mCache(),
  mConnectivityMapPtr(),
  mMeshPtr(new MeshType()) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateBase<Dimension>::
StateBase(const StateBase<Dimension>& rhs):
  mStorage(rhs.mStorage),
  mCache(),
  mNodeListPtrs(rhs.mNodeListPtrs),
  mConnectivityMapPtr(rhs.mConnectivityMapPtr),
  mMeshPtr(rhs.mMeshPtr) {
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateBase<Dimension>::
~StateBase() {
}

//------------------------------------------------------------------------------
// Assignment.
//------------------------------------------------------------------------------
template<typename Dimension>
StateBase<Dimension>&
StateBase<Dimension>::
operator=(const StateBase<Dimension>& rhs) {
  if (this != &rhs) {
    mStorage = rhs.mStorage;
    mCache = CacheType();
    mNodeListPtrs = rhs.mNodeListPtrs;
    mConnectivityMapPtr = rhs.mConnectivityMapPtr;
    mMeshPtr = rhs.mMeshPtr;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Test if the given field name is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
fieldNameRegistered(const FieldName& name) const {
  KeyType fieldName, nodeListName;
  auto itr = mStorage.begin();
  while (itr != mStorage.end()) {
    splitFieldKey(itr->first, fieldName, nodeListName);
    if (fieldName == name) return true;
    ++itr;
  }
  return false;
}

//------------------------------------------------------------------------------
// Test if the internal state is equal.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
operator==(const StateBase<Dimension>& rhs) const {
  if (mStorage.size() != rhs.mStorage.size()) {
    cerr << "Storage sizes don't match." << endl;
    return false;
  }
  vector<KeyType> lhsKeys = keys();
  vector<KeyType> rhsKeys = rhs.keys();
  if (lhsKeys.size() != rhsKeys.size()) {
    cerr << "Keys sizes don't match." << endl;
    return false;
  }
  sort(lhsKeys.begin(), lhsKeys.end());
  sort(rhsKeys.begin(), rhsKeys.end());
  if (lhsKeys != rhsKeys) {
    cerr << "Keys don't match." << endl;
    return false;
  }

  // Walk the keys, and rely on the virtual overloaded 
  // Field::operator==(FieldBase) to do the right thing!
  // We are also relying here on the fact that std::map with a given
  // set of keys will always result in the same order.
  bool result = true;
  typename StorageType::const_iterator lhsItr, rhsItr;
  for (rhsItr = rhs.mStorage.begin(), lhsItr = mStorage.begin();
       rhsItr != rhs.mStorage.end();
       ++rhsItr, ++lhsItr) {
    if (not (*(lhsItr->second) == *(rhsItr->second))) {
      cerr << "Fields for " << lhsItr->first <<  " don't match." << endl;
      result = false;
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Enroll an external mesh.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateBase<Dimension>::
enrollMesh(typename StateBase<Dimension>::MeshPtr meshPtr) {
  mMeshPtr = meshPtr;
}

//------------------------------------------------------------------------------
// Return the full set of known keys.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename StateBase<Dimension>::KeyType>
StateBase<Dimension>::
keys() const {
  vector<KeyType> result;
  result.reserve(mStorage.size());
  for (typename StorageType::const_iterator itr = mStorage.begin();
       itr != mStorage.end();
       ++itr) result.push_back(itr->first);
  ENSURE(result.size() == mStorage.size());
  return result;
}

//------------------------------------------------------------------------------
// Return the set of registered Field names.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename FieldBase<Dimension>::FieldName>
StateBase<Dimension>::
fieldKeys() const {
  KeyType fieldName, nodeListName;
  vector<typename FieldBase<Dimension>::FieldName> result;
  result.reserve(mStorage.size());
  for (typename StorageType::const_iterator itr = mStorage.begin();
       itr != mStorage.end();
       ++itr) {
    splitFieldKey(itr->first, fieldName, nodeListName);
    if (fieldName != "" and nodeListName != "") result.push_back(fieldName);
  }

  // Remove any duplicates.  This will happen when we've stored the same field
  // for different NodeLists.
  sort(result.begin(), result.end());
  result.erase(unique(result.begin(), result.end()), result.end());

  return result;
}

//------------------------------------------------------------------------------
// Enroll a ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateBase<Dimension>::
enrollConnectivityMap(typename StateBase<Dimension>::ConnectivityMapPtr connectivityMapPtr) {
  mConnectivityMapPtr = connectivityMapPtr;
}

//------------------------------------------------------------------------------
// Return the ConnectivityMap.
//------------------------------------------------------------------------------
template<typename Dimension>
const typename StateBase<Dimension>::ConnectivityMapType&
StateBase<Dimension>::
connectivityMap() const {
  REQUIRE(mConnectivityMapPtr.use_count() != 0);
  return *mConnectivityMapPtr;
}

//------------------------------------------------------------------------------
// Return the mesh. (non-const version)
//------------------------------------------------------------------------------
template<typename Dimension>
Mesh<Dimension>&
StateBase<Dimension>::
mesh() {
  REQUIRE(meshRegistered());
  return *mMeshPtr;
}

//------------------------------------------------------------------------------
// Return the mesh. (const version)
//------------------------------------------------------------------------------
template<typename Dimension>
const Mesh<Dimension>&
StateBase<Dimension>::
mesh() const {
  REQUIRE(meshRegistered());
  return *mMeshPtr;
}

//------------------------------------------------------------------------------
// Assign the state data of this State object to be equal to the values
// in another.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateBase<Dimension>::
assign(const StateBase<Dimension>& rhs) {

  // Extract the keys for each state, and verify they line up.
  REQUIRE(mStorage.size() == rhs.mStorage.size());
  vector<KeyType> lhsKeys = keys();
  vector<KeyType> rhsKeys = rhs.keys();
  REQUIRE(lhsKeys.size() == rhsKeys.size());
  sort(lhsKeys.begin(), lhsKeys.end());
  sort(rhsKeys.begin(), rhsKeys.end());
  REQUIRE(lhsKeys == rhsKeys);

  // Walk the keys, and rely on the virtual overloaded 
  // Field::operator=(FieldBase) to do the right thing!
  for (typename StorageType::const_iterator itr = rhs.mStorage.begin();
       itr != rhs.mStorage.end();
       ++itr) {
    *(mStorage[itr->first]) = *(itr->second);
  }

  // Copy the connectivity (by reference).  This thing is too
  // big to carry around separate copies!
  if (rhs.mConnectivityMapPtr != NULL) {
    mConnectivityMapPtr = rhs.mConnectivityMapPtr;
  } else {
    mConnectivityMapPtr = ConnectivityMapPtr();
  }

  // Copy the mesh.
  if (rhs.mMeshPtr != NULL) {
    mMeshPtr = MeshPtr(new MeshType());
    *mMeshPtr = *(rhs.mMeshPtr);
  } else {
    mMeshPtr = MeshPtr();
  }
}

//------------------------------------------------------------------------------
// Force the state fields to be copied to local storage.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateBase<Dimension>::
copyState() {

  // Remove any pre-existing stuff.
  mCache = CacheType();

  // Walk the registered state and copy it to our local cache.
  for (typename StorageType::iterator itr = mStorage.begin();
       itr != mStorage.end();
       ++itr) {
    mCache.push_back(itr->second->clone());
    itr->second = mCache.back().get();
  }
}

}

