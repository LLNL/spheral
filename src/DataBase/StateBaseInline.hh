#include "boost/algorithm/string.hpp"
#include "NodeList/NodeListRegistrar.hh"
#include "Mesh/generateMesh.hh"
#include "Mesh/Mesh.hh"
#include "Utilities/DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Construct the lookup key for the given field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename StateBase<Dimension>::KeyType
StateBase<Dimension>::
key(const FieldSpace::FieldBase<Dimension>& field) {
  return buildFieldKey(field.name(), field.nodeListPtr()->name());
}

//------------------------------------------------------------------------------
// Test if the given key is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
StateBase<Dimension>::
registered(const StateBase<Dimension>::KeyType& key) const {
  return (mStorage.find(key) != mStorage.end());
}

//------------------------------------------------------------------------------
// Test if the given field is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
StateBase<Dimension>::
registered(const FieldSpace::FieldBase<Dimension>& field) const {
  const KeyType key = this->key(field);
  typename StorageType::const_iterator itr = mStorage.find(key);
  return (itr != mStorage.end());
}

//------------------------------------------------------------------------------
// Test if the given field name is registered.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
StateBase<Dimension>::
fieldNameRegistered(const FieldName& name) const {
  KeyType fieldName, nodeListName;
  typename StorageType::const_iterator itr = mStorage.begin();
  while (itr != mStorage.end()) {
    splitFieldKey(itr->first, fieldName, nodeListName);
    if (fieldName == fieldName) return true;
    ++itr;
  }
  return false;
}

//------------------------------------------------------------------------------
// Add a field.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
StateBase<Dimension>::
enroll(FieldSpace::FieldBase<Dimension>& field) {
  const KeyType key = this->key(field);
  mStorage[key] = &field;
  mNodeListPtrs.insert(field.nodeListPtr());
//   std::cerr << "StateBase::enroll storing field:  " << key << " at " << &field << std::endl;
  ENSURE(find(mNodeListPtrs.begin(), mNodeListPtrs.end(), field.nodeListPtr()) != mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Add the fields from a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
void
StateBase<Dimension>::
enrollFieldList(FieldSpace::FieldList<Dimension, Value>& fieldList) {
  for (typename FieldSpace::FieldList<Dimension, Value>::const_iterator itr = fieldList.begin();
       itr != fieldList.end();
       ++itr) {
    this->enroll(**itr);
  }
}

//------------------------------------------------------------------------------
// Return the Field for the given key.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
FieldSpace::Field<Dimension, Value>&
StateBase<Dimension>::
field(const typename StateBase<Dimension>::KeyType& key, 
      const Value& dummy) const {
  try {
    FieldSpace::Field<Dimension, Value>& result = *dynamic_cast<FieldSpace::Field<Dimension, Value>*>(mStorage.find(key)->second);
    return result;
  } catch (const std::bad_cast&) {
    VERIFY2(false, "StateBase::field ERROR: bad field value type.");
  }
}

//------------------------------------------------------------------------------
// Return a FieldList containing all registered fields of the given name.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
FieldSpace::FieldList<Dimension, Value>
StateBase<Dimension>::
fields(const std::string& name, const Value& dummy) const {
  FieldSpace::FieldList<Dimension, Value> result;
  KeyType fieldName, nodeListName;
  for (typename StorageType::const_iterator itr = mStorage.begin();
       itr != mStorage.end();
       ++itr) {
    splitFieldKey(itr->first, fieldName, nodeListName);
    if (fieldName == name) {
      CHECK(nodeListName != "");
      result.appendField(this->field<Value>(itr->first, dummy));
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Return all the Fields of the given Value element type.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename Value>
inline
std::vector<FieldSpace::Field<Dimension, Value>*>
StateBase<Dimension>::
allFields(const Value& dummy) const {
  std::vector<FieldSpace::Field<Dimension, Value>*> result;
  KeyType fieldName, nodeListName;
  for (typename StorageType::const_iterator itr = mStorage.begin();
       itr != mStorage.end();
       ++itr) {
    try {
      FieldSpace::Field<Dimension, Value>* ptr = dynamic_cast<FieldSpace::Field<Dimension, Value>*>(itr->second);
      if (ptr != 0) result.push_back(ptr);
    } catch (const std::bad_cast&) {
      // The field must have been the wrong type.
    }
  }
  return result;
}

//------------------------------------------------------------------------------
// Create a mesh based on the known NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename BoundaryIterator>
inline
void
StateBase<Dimension>::
generateMesh(const Vector& xmin,
             const Vector& xmax,
             const bool generateVoid,
             const bool generateParallelConnectivity,
             const double voidThreshold,
             const BoundaryIterator boundaryBegin,
             const BoundaryIterator boundaryEnd) {
  using NodeSpace::NodeList;
  using std::vector;
  using std::sort;
  NodeList<Dimension> voidNodes("void", 0, 0);
  vector<const NodeList<Dimension>*> nodeLists(mNodeListPtrs.begin(), mNodeListPtrs.end());
  nodeLists.push_back(&voidNodes);
  sort(nodeLists.begin(), nodeLists.end(), typename NodeListRegistrar<Dimension>::NodeListComparator());
  mMeshPtr = boost::shared_ptr<MeshType>(new MeshType());
  MeshSpace::generateMesh<Dimension, 
                          typename vector<const NodeList<Dimension>*>::iterator,
                          BoundaryIterator>
    (nodeLists.begin(), nodeLists.end(),
     boundaryBegin, boundaryEnd,
     xmin, xmax, 
     generateVoid,
     generateParallelConnectivity,
     (generateVoid ? false : true),                           // remove boundary zones
     voidThreshold,
     *mMeshPtr,
     voidNodes);
}

//------------------------------------------------------------------------------
// Test if a mesh is currently available.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
bool
StateBase<Dimension>::
meshRegistered() const {
  return (mMeshPtr.use_count() != 0);
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
