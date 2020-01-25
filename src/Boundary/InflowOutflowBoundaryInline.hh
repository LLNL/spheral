#include <typeinfo>

namespace Spheral {

//------------------------------------------------------------------------------
// The effective timestep vote
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
InflowOutflowBoundary<Dimension>::
dtmin() const {
  return mDT;
}

//------------------------------------------------------------------------------
// Return the DataBase this Boundary is defined on.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const DataBase<Dimension>&
InflowOutflowBoundary<Dimension>::
dataBase() const {
  return mDataBase;
}

//------------------------------------------------------------------------------
// The plane for points to enter by.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const GeomPlane<Dimension>& 
InflowOutflowBoundary<Dimension>::
plane() const {
  return mPlane;
}

//------------------------------------------------------------------------------
// Return the number of nodes this Boundary is going to create.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
InflowOutflowBoundary<Dimension>::
numInflowNodes(const NodeList<Dimension>& nodeList) const {
  const auto itr = mNumInflowNodes.find(nodeList.name());
  VERIFY2(itr != mNumInflowNodes.end(), "InflowOutflowBoundary::numInflowNodes no entry for " << nodeList.name());
  return itr->second;
}

//------------------------------------------------------------------------------
// Access the stored template field values for ghost points.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
std::vector<DataType>&
InflowOutflowBoundary<Dimension>::
storedValues(const KeyType key, const DataType& dummy) {
  auto& storage = storageForType(dummy);
  if (storage.find(key) == storage.end()) {
    VERIFY2(false, "InflowOutflowBoundary ERROR: attempt to extract stored value for " << key << ", which is not stored.");
  }
  return storage[key];
}

template<typename Dimension>
template<typename DataType>
inline
std::vector<DataType>&
InflowOutflowBoundary<Dimension>::
storedValues(const Field<Dimension, DataType>& field) {
  const auto key = StateBase<Dimension>::key(field);
  return storedValues(key, DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Set the stored template field values for ghost points.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
InflowOutflowBoundary<Dimension>::
setStoredValues(const KeyType key, const DataType& value) {
  auto& vals = storedValues(key, value);
  for (auto& x: vals) x = value;
}

template<typename Dimension>
template<typename DataType>
inline
void
InflowOutflowBoundary<Dimension>::
setStoredValues(const Field<Dimension, DataType>& field, const DataType& value) {
  auto& vals = storedValues(field);
  for (auto& x: vals) x = value;
}

}
