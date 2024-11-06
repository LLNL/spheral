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

namespace Spheral {

namespace {

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
// Template to downselect comparison in our variant types
//------------------------------------------------------------------------------
template<typename T1>              bool safeCompare(T1& x, const T1& y) { return x == y; }
template<typename T1, typename T2> bool safeCompare(T1& x, const T2& y) { VERIFY2(false, "Bad comparison!"); return false; }

//------------------------------------------------------------------------------
// Template to downselect assignment in our variant types
//------------------------------------------------------------------------------
template<typename T1>              void safeAssign(T1& x, const T1& y) { x = y; }
template<typename T1, typename T2> void safeAssign(T1& x, const T2& y) { VERIFY2(false, "Bad assignment!"); }

// Helper with overloading in std::visit
template<class... Ts> struct overload : Ts... { using Ts::operator()...; };
template<class... Ts> overload(Ts...) -> overload<Ts...>;

}

//------------------------------------------------------------------------------
// Default constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateBase<Dimension>::
StateBase():
  mFieldStorage(),
  mFieldCache(),
  mMiscStorage(),
  mMiscCache(),
  mNodeListPtrs(),
  mConnectivityMapPtr(),
  mMeshPtr(new MeshType()) {
}

//------------------------------------------------------------------------------
// Copy constructor.
//------------------------------------------------------------------------------
template<typename Dimension>
StateBase<Dimension>::
StateBase(const StateBase<Dimension>& rhs):
  mFieldStorage(rhs.mFieldStorage),
  mFieldCache(),
  mMiscStorage(rhs.mMiscStorage),
  mMiscCache(),
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
    mFieldStorage = rhs.mFieldStorage;
    mFieldCache = FieldCacheType();
    mMiscStorage = rhs.mMiscStorage;
    mMiscCache = MiscCacheType();
    mNodeListPtrs = rhs.mNodeListPtrs;
    mConnectivityMapPtr = rhs.mConnectivityMapPtr;
    mMeshPtr = rhs.mMeshPtr;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Test if the internal state is equal.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
operator==(const StateBase<Dimension>& rhs) const {

  // Compare raw sizes
  if (mFieldStorage.size() != rhs.mFieldStorage.size()) {
    cerr << "Field storage sizes don't match." << endl;
    return false;
  }
  if (mMiscStorage.size() != rhs.mMiscStorage.size()) {
    cerr << "Miscellaneous storage sizes don't match." << endl;
    return false;
  }

  // Keys
  auto lhsKeys = keys();
  auto rhsKeys = rhs.keys();
  if (lhsKeys != rhsKeys) {
    cerr << "Keys don't match." << endl;
    return false;
  }

  // Compare fields
  {
    auto lhsitr = mFieldStorage.begin();
    auto rhsitr = rhs.mFieldStorage.begin();
    for (; lhsitr != mFieldStorage.end(); ++lhsitr, ++rhsitr) {
      CHECK(rhsitr != rhs.mFieldStorage.end());
      CHECK(lhsitr->first == rhsitr->first);
      if (*(lhsitr->second) != *(rhsitr->second)) {
        cerr << "Fields don't match for key " << lhsitr->first << endl;
        return false;
      }
    }
  }

  // Compare the miscellaneous objects
  {
    auto lhsitr = mMiscStorage.begin();
    auto rhsitr = rhs.mMiscStorage.begin();
    for (; lhsitr != mMiscStorage.end(); ++lhsitr, ++rhsitr) {
      CHECK(rhsitr != rhs.mMiscStorage.end());
      CHECK(lhsitr->first == rhsitr->first);
      auto result = std::visit([](auto& x, auto& y) -> bool { return safeCompare(x, y); }, *(lhsitr->second), *(rhsitr->second));
      if (not result) {
        cerr << "State does not match for key " << lhsitr->first << endl;
        return false;
      }
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
  mFieldStorage[key] = &field;
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
  mFieldStorage[key] = fieldPtr.get();
  mNodeListPtrs.insert(fieldPtr->nodeListPtr());
  mFieldCache.push_back(fieldPtr);
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
  return (mFieldStorage.find(key) != mFieldStorage.end() or
          mMiscStorage.find(key) != mMiscStorage.end());
}

//------------------------------------------------------------------------------
// Test if the given field is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
bool
StateBase<Dimension>::
registered(const FieldBase<Dimension>& field) const {
  const auto key = this->key(field);
  return mFieldStorage.find(key) != mFieldStorage.end();
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
  for (auto [key, valptr]: mFieldStorage) {
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
  for (auto itr = mFieldStorage.begin(); itr != mFieldStorage.end(); ++itr) result.push_back(itr->first);
  for (auto itr = mMiscStorage.begin(); itr != mMiscStorage.end(); ++itr) result.push_back(itr->first);
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
  for (auto itr = mFieldStorage.begin(); itr != mFieldStorage.end(); ++itr) result.push_back(itr->first);
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
  for (auto itr = mMiscStorage.begin(); itr != mMiscStorage.end(); ++itr) result.push_back(itr->first);
  return result;
}

//------------------------------------------------------------------------------
// Return the set of registered Field names.
//------------------------------------------------------------------------------
template<typename Dimension>
std::vector<typename FieldBase<Dimension>::FieldName>
StateBase<Dimension>::
fieldNames() const {
  KeyType fieldName, nodeListName;
  vector<FieldName> result;
  for (auto itr = mFieldStorage.begin(); itr != mFieldStorage.end(); ++itr) {
    splitFieldKey(itr->first, fieldName, nodeListName);
    result.push_back(fieldName);
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

  // Fields
  {
    CHECK(mFieldStorage.size() == rhs.mFieldStorage.size());
    auto lhsitr = mFieldStorage.begin();
    auto rhsitr = rhs.mFieldStorage.begin();
    for (; lhsitr != mFieldStorage.end(); ++lhsitr, ++rhsitr) {
      CHECK(rhsitr != rhs.mFieldStorage.end());
      CHECK(lhsitr->first == rhsitr->first);
      *(lhsitr->second) = *(rhsitr->second);
    }
  }

  // Miscellaneous state
  {
    // Depend on assignment working for our AllowedTypes
    CHECK(mMiscStorage.size() == rhs.mMiscStorage.size());
    auto lhsitr = mMiscStorage.begin();
    auto rhsitr = rhs.mMiscStorage.begin();
    for (; lhsitr != mMiscStorage.end(); ++lhsitr, ++rhsitr) {
      CHECK(rhsitr != rhs.mMiscStorage.end());
      CHECK(lhsitr->first == rhsitr->first);
      std::visit([](auto& lhsval, auto& rhsval) { safeAssign(lhsval, rhsval); }, *(lhsitr->second), *(rhsitr->second));
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
  mFieldCache = FieldCacheType();
  mMiscCache = MiscCacheType();

  // Fields
  for (auto itr = mFieldStorage.begin(); itr != mFieldStorage.end(); ++itr) {
    auto clone = itr->second->clone();
    mFieldCache.push_back(clone);
    itr->second = clone.get();
  }

  // Misc
  for (auto itr = mMiscStorage.begin(); itr != mMiscStorage.end(); ++itr) {
    std::visit(overload{[](const Scalar& x)                       { return std::make_shared<AllowedType>(x); },
                        [](const Vector& x)                       { return std::make_shared<AllowedType>(x); },
                        [](const Tensor& x)                       { return std::make_shared<AllowedType>(x); },
                        [](const SymTensor& x)                    { return std::make_shared<AllowedType>(x); },
                        [](const vector<Scalar>& x)               { return std::make_shared<AllowedType>(x); },
                        [](const vector<Vector>& x)               { return std::make_shared<AllowedType>(x); },
                        [](const vector<Tensor>& x)               { return std::make_shared<AllowedType>(x); },
                        [](const vector<SymTensor>& x)            { return std::make_shared<AllowedType>(x); },
                        [](const set<int>& x)                     { return std::make_shared<AllowedType>(x); },
                        [](const set<RKOrder>& x)                 { return std::make_shared<AllowedType>(x); },
                        [](const ReproducingKernel<Dimension>& x) { return std::make_shared<AllowedType>(x); }
      }, *(itr->second));
    // [&](auto* xptr) {
    //              // auto clone = makeClone(*xptr);
    //              auto clone = std::shared_ptr<AllowedType>(makeClone(*xptr));  // new AllowedType(*xptr));
    //              // auto clone = std::make_shared<AllowedType>(*xptr);
    //              // mMiscCache.push_back(clone);
    //              // itr->second = clone.get();
    //            }, itr->second);
  }
}

}

