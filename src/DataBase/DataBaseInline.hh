#include "Utilities/DBC.hh"

namespace Spheral {
namespace DataBaseSpace {

//------------------------------------------------------------------------------
// Number of NodeLists registered with the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
DataBase<Dimension>::numNodeLists() const {
  return mNodeListPtrs.size();
}

//------------------------------------------------------------------------------
// Number of FluidNodeLists registered with the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
DataBase<Dimension>::numFluidNodeLists() const {
  return mFluidNodeListPtrs.size();
}

//------------------------------------------------------------------------------
// Numbers of nodes in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
int
DataBase<Dimension>::numInternalNodes() const {
  int result = 0;
  for (ConstNodeListIterator nodeListItr = nodeListBegin();
       nodeListItr != nodeListEnd();
       ++nodeListItr) result += (*nodeListItr)->numInternalNodes();
  return result;
}

template<typename Dimension>
inline
int
DataBase<Dimension>::numGhostNodes() const {
  int result = 0;
  for (ConstNodeListIterator nodeListItr = nodeListBegin();
       nodeListItr != nodeListEnd();
       ++nodeListItr) result += (*nodeListItr)->numGhostNodes();
  return result;
}

template<typename Dimension>
inline
int
DataBase<Dimension>::numNodes() const {
  int result = 0;
  for (ConstNodeListIterator nodeListItr = nodeListBegin();
       nodeListItr != nodeListEnd();
       ++nodeListItr) result += (*nodeListItr)->numNodes();
  return result;
}

//------------------------------------------------------------------------------
// Standard STL like iterators for NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename DataBase<Dimension>::NodeListIterator
DataBase<Dimension>::nodeListBegin() {
  return mNodeListPtrs.begin();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::NodeListIterator
DataBase<Dimension>::nodeListEnd() {
  return mNodeListPtrs.end();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConstNodeListIterator
DataBase<Dimension>::nodeListBegin() const {
  return mNodeListPtrs.begin();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConstNodeListIterator
DataBase<Dimension>::nodeListEnd() const {
  return mNodeListPtrs.end();
}

//------------------------------------------------------------------------------
// Standard STL like iterators for FluidNodeLists.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename DataBase<Dimension>::FluidNodeListIterator
DataBase<Dimension>::fluidNodeListBegin() {
  return mFluidNodeListPtrs.begin();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::FluidNodeListIterator
DataBase<Dimension>::fluidNodeListEnd() {
  return mFluidNodeListPtrs.end();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConstFluidNodeListIterator
DataBase<Dimension>::fluidNodeListBegin() const {
  return mFluidNodeListPtrs.begin();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConstFluidNodeListIterator
DataBase<Dimension>::fluidNodeListEnd() const {
  return mFluidNodeListPtrs.end();
}

//------------------------------------------------------------------------------
// Standard STL like iterators for FluidNodeLists, but over NodeList
// base type.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename DataBase<Dimension>::NodeListIterator
DataBase<Dimension>::fluidNodeListAsNodeListBegin() {
  return mFluidNodeListAsNodeListPtrs.begin();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::NodeListIterator
DataBase<Dimension>::fluidNodeListAsNodeListEnd() {
  return mFluidNodeListAsNodeListPtrs.end();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConstNodeListIterator
DataBase<Dimension>::fluidNodeListAsNodeListBegin() const {
  return mFluidNodeListAsNodeListPtrs.begin();
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConstNodeListIterator
DataBase<Dimension>::fluidNodeListAsNodeListEnd() const {
  return mFluidNodeListAsNodeListPtrs.end();
}

//------------------------------------------------------------------------------
// Get the current connectivity map.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
const NeighborSpace::ConnectivityMap<Dimension>&
DataBase<Dimension>::
connectivityMap() const {
  VERIFY2(mConnectivityMapPtr.use_count() != 0,
          "DataBase::connectivityMap ERROR -- need to ensure ConnectivityMap is constructed before calling this method.");
  return *mConnectivityMapPtr;
}

template<typename Dimension>
inline
const NeighborSpace::ConnectivityMap<Dimension>&
DataBase<Dimension>::
connectivityMap(const bool buildConnectivityMap) const {
  if (mConnectivityMapPtr.use_count() == 0) this->updateConnectivityMap(buildConnectivityMap);
  return *mConnectivityMapPtr;
}

template<typename Dimension>
inline
typename DataBase<Dimension>::ConnectivityMapPtr
DataBase<Dimension>::
connectivityMapPtr(const bool buildConnectivityMap) const {
  if (mConnectivityMapPtr.use_count() == 0) this->updateConnectivityMap(buildConnectivityMap);
  return mConnectivityMapPtr;
}

//------------------------------------------------------------------------------
// Convenience method to construct a new FieldList with a Field for every
// NodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
FieldSpace::FieldList<Dimension, DataType>
DataBase<Dimension>::
newGlobalFieldList(const DataType value,
                   const typename FieldSpace::Field<Dimension, DataType>::FieldName name) const {
  FieldSpace::FieldList<Dimension, DataType> result(FieldSpace::Copy);
  for (ConstNodeListIterator nodeListItr = nodeListBegin();
       nodeListItr != nodeListEnd();
       ++nodeListItr) {
    result.appendNewField(name, **nodeListItr, value);
  }

  ENSURE(result.numFields() == numNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Convenience method to construct a new FieldList with a Field for every
// FluidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
FieldSpace::FieldList<Dimension, DataType>
DataBase<Dimension>::
newFluidFieldList(const DataType value,
                  const typename FieldSpace::Field<Dimension, DataType>::FieldName name) const {
  FieldSpace::FieldList<Dimension, DataType> result(FieldSpace::Copy);
  for (ConstFluidNodeListIterator nodeListItr = fluidNodeListBegin();
       nodeListItr != fluidNodeListEnd();
       ++nodeListItr) {
    result.appendNewField(name, **nodeListItr, value);
  }

  ENSURE(result.numFields() == numFluidNodeLists());
  return result;
}

//------------------------------------------------------------------------------
// Convenience method to resize a FieldList such that it has a Field for every
// NodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeGlobalFieldList(FieldSpace::FieldList<Dimension, DataType>& fieldList,
                      const DataType value,
                      const typename FieldSpace::Field<Dimension, DataType>::FieldName name,
                      const bool resetValues) const {
  VERIFY((fieldList.storageType() == FieldSpace::Copy));

  // First check if it's necessary to resize the FieldList.
  bool reinitialize = fieldList.numFields() != numNodeLists();
  ConstNodeListIterator nodeListItr = nodeListBegin();
  typename FieldSpace::FieldList<Dimension, DataType>::const_iterator itr = fieldList.begin();
  while (!reinitialize && 
         nodeListItr != nodeListEnd() &&
         itr != fieldList.end()) {
    reinitialize = (*itr)->nodeListPtr() != *nodeListItr;
    ++nodeListItr;
    ++itr;
  }

  if (reinitialize) {
    fieldList = FieldSpace::FieldList<Dimension, DataType>(fieldList.storageType());
    for (ConstNodeListIterator nodeListItr = nodeListBegin();
         nodeListItr != nodeListEnd();
         ++nodeListItr) fieldList.appendNewField(name, **nodeListItr, value);
  } else if (resetValues) {
    fieldList = value;
  }

  ENSURE(fieldList.numFields() == numNodeLists());
}

//------------------------------------------------------------------------------
// Convenience method to resize a FieldList such that it has a Field for every
// FluidNodeList in the DataBase.
//------------------------------------------------------------------------------
template<typename Dimension>
template<typename DataType>
inline
void
DataBase<Dimension>::
resizeFluidFieldList(FieldSpace::FieldList<Dimension, DataType>& fieldList,
                     const DataType value,
                     const typename FieldSpace::Field<Dimension, DataType>::FieldName name,
                     const bool resetValues) const {
  VERIFY((fieldList.storageType() == FieldSpace::Copy));

  // First check if it's necessary to resize the FieldList.
  bool reinitialize = fieldList.numFields() != numFluidNodeLists();
  ConstFluidNodeListIterator nodeListItr = fluidNodeListBegin();
  typename FieldSpace::FieldList<Dimension, DataType>::const_iterator itr = fieldList.begin();
  while (!reinitialize && 
         nodeListItr != fluidNodeListEnd() &&
         itr != fieldList.end()) {
    reinitialize = (*itr)->nodeListPtr() != *nodeListItr;
    ++nodeListItr;
    ++itr;
  }

  if (reinitialize) {
    fieldList = FieldSpace::FieldList<Dimension, DataType>(fieldList.storageType());
    for (ConstFluidNodeListIterator nodeListItr = fluidNodeListBegin();
         nodeListItr != fluidNodeListEnd();
         ++nodeListItr) fieldList.appendNewField(name, **nodeListItr, value);
  } else if (resetValues) {
    fieldList = value;
  }

  ENSURE(fieldList.numFields() == numFluidNodeLists());
}

}
}
