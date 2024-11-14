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
#include "RK/RKCorrectionParams.hh"
#include "RK/ReproducingKernel.hh"
#include "Utilities/AnyVisitor.hh"
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
using std::sort;
using std::shared_ptr;
using std::make_shared;
using std::any;
using std::any_cast;

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Template for generic cloning during copyState
//------------------------------------------------------------------------------
template<typename T>
void
genericClone(std::any& x,
             const std::string& key,
             typename std::map<std::string, std::any>& storage,
             typename std::list<std::any>& cache) {
  auto clone = std::make_shared<T>(std::any_cast<std::reference_wrapper<T>>(x).get());
  cache.push_back(clone);
  storage[key] = std::ref(*clone);
}

}

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateBase<Dimension>::
StateBase():
  mStorage(),
  mCache(),
  mNodeListPtrs(),
  mConnectivityMapPtr(),
  mMeshPtr(new MeshType()) {
}

//------------------------------------------------------------------------------
// Test if the internal state is equal.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
operator==(const StateBase<Dimension>& rhs) const {

  // Compare raw sizes
  if (mStorage.size() != rhs.mStorage.size()) {
    cerr << "Storage sizes don't match." << endl;
    return false;
  }

  // Keys
  auto lhsKeys = keys();
  auto rhsKeys = rhs.keys();
  if (lhsKeys != rhsKeys) {
    cerr << "Keys don't match." << endl;
    return false;
  }

  // Build up a visitor to compare each type of state data we support holding
  AnyVisitor<bool, const std::any&, const std::any&> EQUAL;
  EQUAL.addVisitor<std::reference_wrapper<FieldBase<Dimension>>>        ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(x).get()         == std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<Scalar>>                      ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<Scalar>>(x).get()                       == std::any_cast<std::reference_wrapper<Scalar>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<Vector>>                      ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<Vector>>(x).get()                       == std::any_cast<std::reference_wrapper<Vector>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<Tensor>>                      ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<Tensor>>(x).get()                       == std::any_cast<std::reference_wrapper<Tensor>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<SymTensor>>                   ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<SymTensor>>(x).get()                    == std::any_cast<std::reference_wrapper<SymTensor>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<vector<Scalar>>>              ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<vector<Scalar>>>(x).get()               == std::any_cast<std::reference_wrapper<vector<Scalar>>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<vector<Vector>>>              ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<vector<Vector>>>(x).get()               == std::any_cast<std::reference_wrapper<vector<Vector>>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<vector<Tensor>>>              ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<vector<Tensor>>>(x).get()               == std::any_cast<std::reference_wrapper<vector<Tensor>>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<vector<SymTensor>>>           ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<vector<SymTensor>>>(x).get()            == std::any_cast<std::reference_wrapper<vector<SymTensor>>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<set<int>>>                    ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<set<int>>>(x).get()                     == std::any_cast<std::reference_wrapper<set<int>>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<set<RKOrder>>>                ([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<set<RKOrder>>>(x).get()                 == std::any_cast<std::reference_wrapper<set<RKOrder>>>(y).get(); });
  EQUAL.addVisitor<std::reference_wrapper<ReproducingKernel<Dimension>>>([](const std::any& x, const std::any& y) -> bool { return std::any_cast<std::reference_wrapper<ReproducingKernel<Dimension>>>(x).get() == std::any_cast<std::reference_wrapper<ReproducingKernel<Dimension>>>(y).get(); });
  
  // Apply the equality visitor to all the stored State data
  auto lhsitr = mStorage.begin();
  auto rhsitr = rhs.mStorage.begin();
  for (; lhsitr != mStorage.end(); ++lhsitr, ++rhsitr) {
    CHECK(rhsitr != rhs.mStorage.end());
    CHECK(lhsitr->first == rhsitr->first);
    if (not EQUAL.visit(lhsitr->second, rhsitr->second)) {
      cerr << "States don't match for key " << lhsitr->first << endl;
      return false;
    }
  }

  return true;
}

//------------------------------------------------------------------------------
// Enroll a Field
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateBase<Dimension>::
enroll(FieldBase<Dimension>& field) {
  const auto key = this->key(field);
  mStorage[key] = std::ref(field);
  mNodeListPtrs.insert(field.nodeListPtr());
  // std::cerr << "StateBase::enroll field:  " << key << " at " << &field << std::endl;
  ENSURE(std::find(mNodeListPtrs.begin(), mNodeListPtrs.end(), field.nodeListPtr()) != mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Enroll a Field (shared_ptr).
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateBase<Dimension>::
enroll(std::shared_ptr<FieldBase<Dimension>>& fieldPtr) {
  const auto key = this->key(*fieldPtr);
  mStorage[key] = std::ref(*fieldPtr);
  mNodeListPtrs.insert(fieldPtr->nodeListPtr());
  mCache.push_back(fieldPtr);
  ENSURE(std::find(mNodeListPtrs.begin(), mNodeListPtrs.end(), fieldPtr->nodeListPtr()) != mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Add the fields from a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
void
StateBase<Dimension>::
enroll(FieldListBase<Dimension>& fieldList) {
  for (auto* fptr: range(fieldList.begin_base(), fieldList.end_base())) {
    this->enroll(*fptr);
  }
}

//------------------------------------------------------------------------------
// Test if the given key is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
registered(const StateBase<Dimension>::KeyType& key) const {
  return mStorage.find(key) != mStorage.end();
}

//------------------------------------------------------------------------------
// Test if the given field is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
registered(const FieldBase<Dimension>& field) const {
  const auto key = this->key(field);
  return this->registered(key);
}

//------------------------------------------------------------------------------
// Test if the given FieldList is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
registered(const FieldListBase<Dimension>& fieldList) const {
  REQUIRE(fieldList.begin_base() != fieldList.end_base());
  return this->registered(**fieldList.begin_base());
}

//------------------------------------------------------------------------------
// Test if the given field name is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
fieldNameRegistered(const FieldName& name) const {
  KeyType fieldName, nodeListName;
  for (auto [key, valptr]: mStorage) {
    splitFieldKey(key, fieldName, nodeListName);
    if (fieldName == name) return true;
  }
  return false;
}

//------------------------------------------------------------------------------
// Return the full set of known keys.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename StateBase<Dimension>::KeyType>
StateBase<Dimension>::
keys() const {
  vector<KeyType> result;
  for (auto itr = mStorage.begin(); itr != mStorage.end(); ++itr) result.push_back(itr->first);
  return result;
}

//------------------------------------------------------------------------------
// Return the full set of Field Keys (mangled with NodeList names)
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename StateBase<Dimension>::KeyType>
StateBase<Dimension>::
fullFieldKeys() const {
  vector<KeyType> result;
  for (auto [key, aref]: mStorage) {
    if (aref.type() == typeid(std::reference_wrapper<FieldBase<Dimension>>)) {
      result.push_back(key);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the set of non-field keys.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename StateBase<Dimension>::KeyType>
StateBase<Dimension>::
miscKeys() const {
  vector<KeyType> result;
  for (auto [key, aref]: mStorage) {
    if (aref.type() != typeid(std::reference_wrapper<FieldBase<Dimension>>)) {
      result.push_back(key);
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Return the set of registered Field names.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename FieldBase<Dimension>::FieldName>
StateBase<Dimension>::
fieldNames() const {
  vector<FieldName> result;
  KeyType fieldName, nodeListName;
  for (auto [key, aref]: mStorage) {
    if (aref.type() == typeid(std::reference_wrapper<FieldBase<Dimension>>)) {
      auto fref = std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(aref);
      splitFieldKey(fref.get().name(), fieldName, nodeListName);
      result.push_back(fieldName);
    }
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
// Return the ConnectivityMap (non-const version)
//------------------------------------------------------------------------------
template<typename Dimension>
typename StateBase<Dimension>::ConnectivityMapType&
StateBase<Dimension>::
connectivityMap() {
  REQUIRE(mConnectivityMapPtr.use_count() != 0);
  return *mConnectivityMapPtr;
}

//------------------------------------------------------------------------------
// Return the ConnectivityMap (const version)
//------------------------------------------------------------------------------
template<typename Dimension>
const typename StateBase<Dimension>::ConnectivityMapType&
StateBase<Dimension>::
connectivityMap() const {
  REQUIRE(mConnectivityMapPtr.use_count() != 0);
  return *mConnectivityMapPtr;
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
// Test if a mesh is currently available.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
meshRegistered() const {
  return (mMeshPtr.use_count() != 0);
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

  // Build a visitor that knows how to assign each of our datatypes
  AnyVisitor<void, std::any&, const std::any&> ASSIGN;
  ASSIGN.addVisitor<std::reference_wrapper<FieldBase<Dimension>>>        ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(x).get()         = std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<Scalar>>                      ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<Scalar>>(x).get()                       = std::any_cast<std::reference_wrapper<Scalar>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<Vector>>                      ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<Vector>>(x).get()                       = std::any_cast<std::reference_wrapper<Vector>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<Tensor>>                      ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<Tensor>>(x).get()                       = std::any_cast<std::reference_wrapper<Tensor>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<SymTensor>>                   ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<SymTensor>>(x).get()                    = std::any_cast<std::reference_wrapper<SymTensor>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<vector<Scalar>>>              ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<vector<Scalar>>>(x).get()               = std::any_cast<std::reference_wrapper<vector<Scalar>>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<vector<Vector>>>              ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<vector<Vector>>>(x).get()               = std::any_cast<std::reference_wrapper<vector<Vector>>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<vector<Tensor>>>              ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<vector<Tensor>>>(x).get()               = std::any_cast<std::reference_wrapper<vector<Tensor>>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<vector<SymTensor>>>           ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<vector<SymTensor>>>(x).get()            = std::any_cast<std::reference_wrapper<vector<SymTensor>>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<set<int>>>                    ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<set<int>>>(x).get()                     = std::any_cast<std::reference_wrapper<set<int>>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<set<RKOrder>>>                ([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<set<RKOrder>>>(x).get()                 = std::any_cast<std::reference_wrapper<set<RKOrder>>>(y).get(); });
  ASSIGN.addVisitor<std::reference_wrapper<ReproducingKernel<Dimension>>>([](std::any& x, const std::any& y) { std::any_cast<std::reference_wrapper<ReproducingKernel<Dimension>>>(x).get() = std::any_cast<std::reference_wrapper<ReproducingKernel<Dimension>>>(y).get(); });

  // Apply the assignment visitor to all the stored State data
  auto lhsitr = mStorage.begin();
  auto rhsitr = rhs.mStorage.begin();
  for (; lhsitr != mStorage.end(); ++lhsitr, ++rhsitr) {
    CHECK(rhsitr != rhs.mStorage.end());
    CHECK(lhsitr->first == rhsitr->first);
    ASSIGN.visit(lhsitr->second, rhsitr->second);
  }

  // Copy the connectivity (by reference).  This thing is too
  // big to carry around separate copies!
  if (rhs.mConnectivityMapPtr != nullptr) {
    mConnectivityMapPtr = rhs.mConnectivityMapPtr;
  } else {
    mConnectivityMapPtr = ConnectivityMapPtr();
  }

  // Copy the mesh.
  if (rhs.mMeshPtr != nullptr) {
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

  // Build a visitor to clone each type of state data
  AnyVisitor<void, std::any&, const KeyType&, StorageType&, CacheType&> CLONE;
  CLONE.addVisitor<std::reference_wrapper<FieldBase<Dimension>>>            ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) {
                                                                               auto clone = std::any_cast<std::reference_wrapper<FieldBase<Dimension>>>(x).get().clone();
                                                                               cache.push_back(clone);
                                                                               storage[key] = std::ref(*clone);
                                                                             });
  CLONE.addVisitor<std::reference_wrapper<Scalar>>                          ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<Scalar>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<Vector>>                          ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<Vector>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<Tensor>>                          ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<Tensor>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<SymTensor>>                       ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<SymTensor>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<vector<Scalar>>>                  ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<Scalar>>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<vector<Vector>>>                  ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<Vector>>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<vector<Tensor>>>                  ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<Tensor>>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<vector<SymTensor>>>               ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<SymTensor>>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<set<int>>>                        ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<set<int>>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<set<RKOrder>>>                    ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<set<RKOrder>>(x, key, storage, cache); });
  CLONE.addVisitor<std::reference_wrapper<ReproducingKernel<Dimension>>>    ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<ReproducingKernel<Dimension>>(x, key, storage, cache); });

  // Clone all our stored data to cache
  for (auto& [key, anyval]: mStorage) {
    CLONE.visit(anyval, key, mStorage, mCache);
  }
}

}

