#include <typeinfo>

namespace Spheral {

//------------------------------------------------------------------------------
// Return the number of nodes this Boundary is going to create.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
InflowBoundary<Dimension>::
numInflowNodes() const {
  return mNumInflowNodes;
}

//------------------------------------------------------------------------------
// Return the NodeList this Boundary is defined on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NodeList<Dimension>&
InflowBoundary<Dimension>::
nodeList() const {
  CHECK(mNodeListPtr != 0);
  return *mNodeListPtr;
}

//------------------------------------------------------------------------------
// The plane for points to enter by.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>& 
InflowBoundary<Dimension>::
plane() const {
  return mPlane;
}

//------------------------------------------------------------------------------
// Access the stored template field values for ghost points.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
std::vector<DataType>&
InflowBoundary<Dimension>::
storedValues(const std::string fieldName, const DataType& dummy) {
  auto& storage = storageForType(dummy);
  const auto key = StateBase<Dimension>::buildFieldKey(fieldName, mNodeListPtr->name());
  if (storage.find(key) == storage.end()) {
    VERIFY2(false, "InflowBoundary ERROR: attempt to extract stored value for " << fieldName << ", which is not stored.");
  }
  return storage[key];
}

template<typename Dimension>
template<typename DataType>
inline
std::vector<DataType>&
InflowBoundary<Dimension>::
storedValues(const Field<Dimension, DataType>& field) {
  return storedValues(field.name(), DataTypeTraits<DataType>::zero());
}

}
