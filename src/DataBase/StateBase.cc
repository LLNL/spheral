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
// Collect visitor methods to apply to std::any object holders
//------------------------------------------------------------------------------
// 2 args
template<typename RETURNT, typename ARG1, typename ARG2>
class AnyVisitor2 {
public:
  using VisitorFunc = std::function<RETURNT (ARG1, ARG2)>;

  RETURNT visit(ARG1 value1, ARG2 value2) const {
    auto it = mVisitors.find(std::type_index(value1.type()));
    if (it != mVisitors.end()) {
      return it->second(value1, value2);
    }
    VERIFY2(false, "AnyVisitor ERROR in StateBase: unable to process unknown data");
  }

  template<typename T>
  void addVisitor(VisitorFunc visitor) {
    mVisitors[std::type_index(typeid(T))] = visitor;
  }


private:
  std::unordered_map<std::type_index, VisitorFunc> mVisitors;
};

//..............................................................................
// 4 args
template<typename RETURNT, typename ARG1, typename ARG2, typename ARG3, typename ARG4>
class AnyVisitor4 {
public:
  using VisitorFunc = std::function<RETURNT (ARG1, ARG2, ARG3, ARG4)>;

  RETURNT visit(ARG1 value1, ARG2 value2, ARG3 value3, ARG4 value4) const {
    auto it = mVisitors.find(std::type_index(value1.type()));
    if (it != mVisitors.end()) {
      return it->second(value1, value2, value3, value4);
    }
    VERIFY2(false, "AnyVisitor ERROR in StateBase: unable to process unknown data");
  }

  template<typename T>
  void addVisitor(VisitorFunc visitor) {
    mVisitors[std::type_index(typeid(T))] = visitor;
  }


private:
  std::unordered_map<std::type_index, VisitorFunc> mVisitors;
};

// //------------------------------------------------------------------------------
// // Helper for copying a type, used in copyState
// //------------------------------------------------------------------------------
// template<typename T>
// T*
// extractType(boost::any& anyT) {
//   try {
//     T* result = boost::any_cast<T*>(anyT);
//     return result;
//   } catch (boost::any_cast_error) {
//     return NULL;
//   }
// }

//------------------------------------------------------------------------------
// Template for generic cloning during copyState
//------------------------------------------------------------------------------
template<typename T>
void
genericClone(std::any& x,
             const std::string& key,
             typename std::map<std::string, std::any>& storage,
             typename std::list<std::any>& cache) {
  auto clone = std::make_shared<T>(*std::any_cast<T*>(x));
  cache.push_back(clone);
  storage[key] = clone.get();
}

//------------------------------------------------------------------------------
// Template to downselect comparison in our variant types
//------------------------------------------------------------------------------
template<typename T1>              bool safeCompare(T1& x, const T1& y) { return x == y; }
template<typename T1, typename T2> bool safeCompare(T1& x, const T2& y) { VERIFY2(false, "Bad comparison!"); return false; }

//------------------------------------------------------------------------------
// Template to downselect assignment in our variant types
//------------------------------------------------------------------------------
template<typename T1>              void safeAssign(T1& x, const T1& y) { x = y; }
template<typename T1, typename T2> void safeAssign(T1& x, const T2& y) { VERIFY2(false, "Bad assignment!"); }

template<typename T1>              T1& safePointer(T1* xptr, const T1* yptr) { return yptr; }
template<typename T1, typename T2> T1& safePointer(T1* xptr, const T2* yptr) { VERIFY2(false, "Bad assignment!"); return xptr; }

//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template<typename ResultT, typename T1> std::shared_ptr<ResultT> safeClone(const T1& x, const ResultT& dummy) { return std::make_shared<ResultT>(x); }

//------------------------------------------------------------------------------
// Helper with overloading in std::visit
//------------------------------------------------------------------------------
template<class... Ts> struct overload : Ts... { using Ts::operator()...; };
template<class... Ts> overload(Ts...) -> overload<Ts...>;

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
  AnyVisitor2<bool, const std::any&, const std::any&> EQUAL;
  EQUAL.addVisitor<FieldBase<Dimension>*>        ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<FieldBase<Dimension>*>(x)         == *std::any_cast<FieldBase<Dimension>*>(y); });
  EQUAL.addVisitor<Scalar*>                      ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<Scalar*>(x)                       == *std::any_cast<Scalar*>(y); });
  EQUAL.addVisitor<Vector*>                      ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<Vector*>(x)                       == *std::any_cast<Vector*>(y); });
  EQUAL.addVisitor<Tensor*>                      ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<Tensor*>(x)                       == *std::any_cast<Tensor*>(y); });
  EQUAL.addVisitor<SymTensor*>                   ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<SymTensor*>(x)                    == *std::any_cast<SymTensor*>(y); });
  EQUAL.addVisitor<vector<Scalar>*>              ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<vector<Scalar>*>(x)               == *std::any_cast<vector<Scalar>*>(y); });
  EQUAL.addVisitor<vector<Vector>*>              ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<vector<Vector>*>(x)               == *std::any_cast<vector<Vector>*>(y); });
  EQUAL.addVisitor<vector<Tensor>*>              ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<vector<Tensor>*>(x)               == *std::any_cast<vector<Tensor>*>(y); });
  EQUAL.addVisitor<vector<SymTensor>*>           ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<vector<SymTensor>*>(x)            == *std::any_cast<vector<SymTensor>*>(y); });
  EQUAL.addVisitor<set<int>*>                    ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<set<int>*>(x)                     == *std::any_cast<set<int>*>(y); });
  EQUAL.addVisitor<set<RKOrder>*>                ([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<set<RKOrder>*>(x)                 == *std::any_cast<set<RKOrder>*>(y); });
  EQUAL.addVisitor<ReproducingKernel<Dimension>*>([](const std::any& x, const std::any& y) -> bool { return *std::any_cast<ReproducingKernel<Dimension>*>(x) == *std::any_cast<ReproducingKernel<Dimension>*>(y); });
  
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
  mStorage[key] = &field;
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
  mStorage[key] = fieldPtr.get();
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
  for (auto [key, aptr]: mStorage) {
    if (std::any_cast<FieldBase<Dimension>*>(aptr) != nullptr) result.push_back(key);
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
  for (auto [key, aptr]: mStorage) {
    if (std::any_cast<FieldBase<Dimension>*>(aptr) == nullptr) result.push_back(key);
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
  for (auto [key, aptr]: mStorage) {
    auto* fptr = std::any_cast<FieldBase<Dimension>*>(aptr);
    if (fptr != nullptr) result.push_back(fptr->name());
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
  AnyVisitor2<void, std::any&, const std::any&> ASSIGN;
  ASSIGN.addVisitor<FieldBase<Dimension>*>        ([](std::any& x, const std::any& y) { *std::any_cast<FieldBase<Dimension>*>(x)         = *std::any_cast<FieldBase<Dimension>*>(y); });
  ASSIGN.addVisitor<Scalar*>                      ([](std::any& x, const std::any& y) { *std::any_cast<Scalar*>(x)                       = *std::any_cast<Scalar*>(y); });
  ASSIGN.addVisitor<Vector*>                      ([](std::any& x, const std::any& y) { *std::any_cast<Vector*>(x)                       = *std::any_cast<Vector*>(y); });
  ASSIGN.addVisitor<Tensor*>                      ([](std::any& x, const std::any& y) { *std::any_cast<Tensor*>(x)                       = *std::any_cast<Tensor*>(y); });
  ASSIGN.addVisitor<SymTensor*>                   ([](std::any& x, const std::any& y) { *std::any_cast<SymTensor*>(x)                    = *std::any_cast<SymTensor*>(y); });
  ASSIGN.addVisitor<vector<Scalar>*>              ([](std::any& x, const std::any& y) { *std::any_cast<vector<Scalar>*>(x)               = *std::any_cast<vector<Scalar>*>(y); });
  ASSIGN.addVisitor<vector<Vector>*>              ([](std::any& x, const std::any& y) { *std::any_cast<vector<Vector>*>(x)               = *std::any_cast<vector<Vector>*>(y); });
  ASSIGN.addVisitor<vector<Tensor>*>              ([](std::any& x, const std::any& y) { *std::any_cast<vector<Tensor>*>(x)               = *std::any_cast<vector<Tensor>*>(y); });
  ASSIGN.addVisitor<vector<SymTensor>*>           ([](std::any& x, const std::any& y) { *std::any_cast<vector<SymTensor>*>(x)            = *std::any_cast<vector<SymTensor>*>(y); });
  ASSIGN.addVisitor<set<int>*>                    ([](std::any& x, const std::any& y) { *std::any_cast<set<int>*>(x)                     = *std::any_cast<set<int>*>(y); });
  ASSIGN.addVisitor<set<RKOrder>*>                ([](std::any& x, const std::any& y) { *std::any_cast<set<RKOrder>*>(x)                 = *std::any_cast<set<RKOrder>*>(y); });
  ASSIGN.addVisitor<ReproducingKernel<Dimension>*>([](std::any& x, const std::any& y) { *std::any_cast<ReproducingKernel<Dimension>*>(x) = *std::any_cast<ReproducingKernel<Dimension>*>(y); });

  // Apply the assignment visitor to all the stored State data
  auto lhsitr = mStorage.begin();
  auto rhsitr = rhs.mStorage.begin();
  for (; lhsitr != mStorage.end(); ++lhsitr, ++rhsitr) {
    CHECK(rhsitr != rhs.mStorage.end());
    CHECK(lhsitr->first == rhsitr->first);
    try {
      ASSIGN.visit(lhsitr->second, rhsitr->second);
    } catch(...) {
      CHECK(false);
    }
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
  AnyVisitor4<void, std::any&, const KeyType&, StorageType&, CacheType&> CLONE;
  CLONE.addVisitor<FieldBase<Dimension>*>            ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) {
                                                        auto clone = std::any_cast<FieldBase<Dimension>*>(x)->clone();
                                                        cache.push_back(clone);
                                                        storage[key] = clone.get();
                                                      });
  CLONE.addVisitor<Scalar*>                          ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<Scalar>(x, key, storage, cache); });
  CLONE.addVisitor<Vector*>                          ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<Vector>(x, key, storage, cache); });
  CLONE.addVisitor<Tensor*>                          ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<Tensor>(x, key, storage, cache); });
  CLONE.addVisitor<SymTensor*>                       ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<SymTensor>(x, key, storage, cache); });
  CLONE.addVisitor<vector<Scalar>*>                  ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<Scalar>>(x, key, storage, cache); });
  CLONE.addVisitor<vector<Vector>*>                  ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<Vector>>(x, key, storage, cache); });
  CLONE.addVisitor<vector<Tensor>*>                  ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<Tensor>>(x, key, storage, cache); });
  CLONE.addVisitor<vector<SymTensor>*>               ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<vector<SymTensor>>(x, key, storage, cache); });
  CLONE.addVisitor<set<int>*>                        ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<set<int>>(x, key, storage, cache); });
  CLONE.addVisitor<set<RKOrder>*>                    ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<set<RKOrder>>(x, key, storage, cache); });
  CLONE.addVisitor<ReproducingKernel<Dimension>*>    ([](std::any& x, const KeyType& key, StorageType& storage, CacheType& cache) { genericClone<ReproducingKernel<Dimension>>(x, key, storage, cache); });

  // Clone all our stored data to cache
  for (auto& [key, anyvalptr]: mStorage) {
    CLONE.visit(anyvalptr, key, mStorage, mCache);
  }
}

}

