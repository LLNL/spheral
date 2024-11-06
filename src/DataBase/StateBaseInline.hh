#include "boost/algorithm/string.hpp"
#include "DataBase/UpdatePolicyBase.hh"
#include "RK/RKCorrectionParams.hh"
#include "RK/ReproducingKernel.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/range.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Enroll an arbitrary type
// Must be one of the supported types in StateBase::AllowedType
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename T>
inline
void
StateBase<Dimension>::
enroll(const KeyType& key, T& thing) {
  mMiscStorage[key] = &thing;
}

//------------------------------------------------------------------------------
// Return the Field for the given key.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
Field<Dimension, Value>&
StateBase<Dimension>::
field(const KeyType& key) const {
  auto itr = mFieldStorage.find(key);
  VERIFY2(itr != mFieldStorage.end(), "StateBase ERROR: failed lookup for Field " << key);
  auto* fbasePtr = itr->second;
  auto* resultPtr = dynamic_cast<Field<Dimension, Value>*>(fbasePtr);
  VERIFY2(resultPtr != nullptr,
          "StateBase::field ERROR: field type incorrect for key " << key);
  return *resultPtr;
}

template<typename Dimension>
template<typename Value>
inline
Field<Dimension, Value>&
StateBase<Dimension>::
field(const typename StateBase<Dimension>::KeyType& key, 
      const Value&) const {
  return this->template field<Value>(key);
}

//------------------------------------------------------------------------------
// Return all the Fields of the given Value element type.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
std::vector<Field<Dimension, Value>*>
StateBase<Dimension>::
allFields(const Value&) const {
  std::vector<Field<Dimension, Value>*> result;
  KeyType fieldName, nodeListName;
  for (auto [key, valptr]: mFieldStorage) {
    auto* ptr = dynamic_cast<Field<Dimension, Value>*>(valptr);
    if (ptr != nullptr) result.push_back(ptr);
  }
  return result;
}

//------------------------------------------------------------------------------
// Return a FieldList containing all registered fields of the given name.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
FieldList<Dimension, Value>
StateBase<Dimension>::
fields(const std::string& name) const {
  FieldList<Dimension, Value> result;
  KeyType fieldName, nodeListName;
  for (auto [key, valptr]: mFieldStorage) {
    splitFieldKey(key, fieldName, nodeListName);
    if (fieldName == name) {
      CHECK(nodeListName != "");
      auto* fptr = dynamic_cast<Field<Dimension, Value>*>(valptr);
      CHECK(valptr != nullptr);
      result.appendField(*fptr);
    }
  }
  return result;
}

template<typename Dimension>
template<typename Value>
inline
FieldList<Dimension, Value>
StateBase<Dimension>::
fields(const std::string& name, const Value& dummy) const {
  return this->template fields<Value>(name);
}

//------------------------------------------------------------------------------
// Extract an arbitrary type
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
Value&
StateBase<Dimension>::
get(const typename StateBase<Dimension>::KeyType& key) const {
  auto itr = mMiscStorage.find(key);
  VERIFY2(itr != mMiscStorage.end(), "StateBase ERROR: failed lookup for key " << key);
  auto* resultPtr = std::get_if<Value>(itr->second);
  VERIFY2(resultPtr != nullptr, "StateBase::get ERROR: unable to extract Value for " << key << "\n");
  return *resultPtr;
}

// Same thing passing a dummy argument to help with template type
template<typename Dimension>
template<typename Value>
inline
Value&
StateBase<Dimension>::
get(const typename StateBase<Dimension>::KeyType& key,
    const Value&) const {
  return this->get<Value>(key);
}

//------------------------------------------------------------------------------
// Assign the Fields matching the given name of this State object to be equal to
// the values in another.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
void
StateBase<Dimension>::
assignFields(const StateBase<Dimension>& rhs, const std::string name) {
  auto lhsfields = this->fields(name, Value());
  auto rhsfields = rhs.fields(name, Value());
  lhsfields.assignFields(rhsfields);
}

//------------------------------------------------------------------------------
// Construct the lookup key for the given field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename StateBase<Dimension>::KeyType
StateBase<Dimension>::
key(const FieldBase<Dimension>& field) {
  return buildFieldKey(field.name(), field.nodeListPtr()->name());
}

//------------------------------------------------------------------------------
// Construct the lookup key for the given FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename StateBase<Dimension>::KeyType
StateBase<Dimension>::
key(const FieldListBase<Dimension>& fieldList) {
  REQUIRE(fieldList.begin_base() != fieldList.end_base());
  return buildFieldKey((*fieldList.begin_base())->name(), UpdatePolicyBase<Dimension>::wildcard());
}

//------------------------------------------------------------------------------
// Internal methods to encode the convention for combining Field and NodeList
// names into a single unique key.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename StateBase<Dimension>::KeyType
StateBase<Dimension>::
buildFieldKey(const std::string& fieldName,
              const std::string& nodeListName) {
  return fieldName + "|" + nodeListName;
}

template<typename Dimension>
inline
void
StateBase<Dimension>::
splitFieldKey(const KeyType& key,
              KeyType& fieldName,
              KeyType& nodeListName) {
  std::vector<std::string> keys;
  boost::split(keys, key, boost::is_any_of("|"));
  if (keys.size() > 1) {
    fieldName = keys[0];
    nodeListName = keys[1];
  } else if (keys.size() == 1) {
    fieldName = keys[0];
    nodeListName = "";
  } else {
    fieldName = "";
    nodeListName = "";
  }
}

}
