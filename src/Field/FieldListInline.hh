// Includes.
#include "Geometry/MathTraits.hh"
#include "NodeIterators.hh"

#include "NodeList/FluidNodeList.hh"
#include "NodeList/NodeListRegistrar.hh"
#include "Neighbor/Neighbor.hh"
#include "Field/Field.hh"
#include "Kernel/TableKernel.hh"

#ifdef USE_MPI
#include "mpi.h"
#include "Utilities/DataTypeTraits.hh"
#include "Utilities/packElement.hh"
#include "Distributed/Communicator.hh"
#endif

#include <algorithm>
#include <limits.h>
#include <float.h>

namespace Spheral {

//------------------------------------------------------------------------------
// Empty constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::FieldList():
  FieldListBase<Dimension>(),
  mFieldPtrs(0),
  mFieldBasePtrs(0),
  mFieldCache(0),
  mStorageType(FieldStorageType::ReferenceFields),
  mNodeListPtrs(0),
  mNodeListIndexMap() {
}

//------------------------------------------------------------------------------
// Construct an empty Field List, with an explict storage type.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::FieldList(FieldStorageType aStorageType):
  FieldListBase<Dimension>(),
  mFieldPtrs(0),
  mFieldBasePtrs(0),
  mFieldCache(0),
  mStorageType(aStorageType),
  mNodeListPtrs(0),
  mNodeListIndexMap() {
}

//------------------------------------------------------------------------------
// Copy constructor.
// Note that the copy constructor copies the data by reference by default.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::
FieldList(const FieldList<Dimension, DataType>& rhs):
  FieldListBase<Dimension>(rhs),
  mFieldPtrs(rhs.mFieldPtrs),
  mFieldBasePtrs(rhs.mFieldBasePtrs),
  mFieldCache(),
  mStorageType(rhs.storageType()),
  mNodeListPtrs(rhs.mNodeListPtrs),
  mNodeListIndexMap(rhs.mNodeListIndexMap) {

  // If we're maintaining Fields by copy, then copy the Field cache.
  if (storageType() == FieldStorageType::CopyFields) {
    CHECK(mFieldCache.size() == 0);
    for (typename FieldCacheType::const_iterator itr = rhs.mFieldCache.begin();
         itr != rhs.mFieldCache.end();
         ++itr) {
      std::shared_ptr<Field<Dimension, DataType> > newField(new Field<Dimension, DataType>(**itr));
      mFieldCache.push_back(newField);
    }

    CHECK(this->size() == mFieldCache.size());
    iterator fieldPtrItr = this->begin();
    typename FieldListBase<Dimension>::iterator fieldBasePtrItr = this->begin_base();
    typename FieldCacheType::iterator fieldCacheItr = mFieldCache.begin();
    for(; fieldPtrItr != this->end(); ++fieldPtrItr, ++fieldBasePtrItr, ++fieldCacheItr) {
      CHECK(fieldCacheItr != mFieldCache.end());
      (*fieldPtrItr) = fieldCacheItr->get();
      (*fieldBasePtrItr) = fieldCacheItr->get();
    }

    CHECK(fieldPtrItr == this->end() &&
          fieldBasePtrItr == this->end_base() &&
          fieldCacheItr == mFieldCache.end());
  }

//   // Register this FieldList with each Field we point to.
//   for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//     registerWithField(**fieldPtrItr);

}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>::~FieldList() {

//   // Unregister this FieldList from each Field we point to.
//   if (storageType() == FieldStorageType::ReferenceFields) {
//     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//       unregisterFromField(**fieldPtrItr);
//   }

}

//------------------------------------------------------------------------------
// Assignment with another FieldList.
// Note that we copy the Field data by Reference.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator=(const FieldList<Dimension, DataType>& rhs) {
  if (this != &rhs) {
    mStorageType = rhs.storageType();
    mNodeListPtrs = rhs.mNodeListPtrs;
    mFieldCache = FieldCacheType();
    mNodeListIndexMap = rhs.mNodeListIndexMap;
    mFieldPtrs = std::vector<ElementType>();
    mFieldBasePtrs = std::vector<BaseElementType>();
    mFieldPtrs.reserve(rhs.size());
    mFieldBasePtrs.reserve(rhs.size());

//     // Unregister from our current set of Fields.
//     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//       unregisterFromField(**fieldPtrItr);

    switch(storageType()) {

    case FieldStorageType::ReferenceFields:
      for (const_iterator fieldPtrItr = rhs.begin(); 
           fieldPtrItr != rhs.end();
           ++fieldPtrItr) {
        mFieldPtrs.push_back(*fieldPtrItr);
        mFieldBasePtrs.push_back(*fieldPtrItr);
      }
      break;

    case FieldStorageType::CopyFields:
      for (typename FieldCacheType::const_iterator itr = rhs.mFieldCache.begin();
         itr != rhs.mFieldCache.end();
         ++itr) {
        std::shared_ptr<Field<Dimension, DataType> > newField(new Field<Dimension, DataType>(**itr));
        mFieldCache.push_back(newField);
      }

      for (typename FieldCacheType::iterator fieldCacheItr = mFieldCache.begin();
           fieldCacheItr != mFieldCache.end();
           ++fieldCacheItr) {
        mFieldPtrs.push_back(fieldCacheItr->get());
        mFieldBasePtrs.push_back(fieldCacheItr->get());
      }
      CHECK(this->size() == mFieldCache.size());
      CHECK(mFieldBasePtrs.size() == mFieldCache.size());
      break;
    }

//     // Register this FieldList with each Field we point to.
//     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//       registerWithField(**fieldPtrItr);

  }
  ENSURE(this->size() == rhs.size());
  return *this;
}

//------------------------------------------------------------------------------
// Assignment with a constant.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator=(const DataType& rhs) {
  for (iterator fieldPtrItr = begin(); fieldPtrItr < end(); ++fieldPtrItr) {
    (*fieldPtrItr)->operator=(rhs);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Return the storage type.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldStorageType
FieldList<Dimension, DataType>::storageType() const {
  return mStorageType;
}

//------------------------------------------------------------------------------
// Force the FieldList to store the fields it points to.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::copyFields() {
  if (storageType() != FieldStorageType::CopyFields) {
    mStorageType = FieldStorageType::CopyFields;

//     // Unregister this FieldList from each Field we currently point to.
//     for (iterator fieldPtrItr = begin(); fieldPtrItr != end(); ++fieldPtrItr) 
//       unregisterFromField(**fieldPtrItr);

    // Store local copies of the fields we're pointing at.
    mFieldCache = FieldCacheType();
    iterator itr = begin();
    typename FieldListBase<Dimension>::iterator baseItr = begin_base();
    for (; itr != end(); ++itr, ++baseItr) {
      std::shared_ptr<Field<Dimension, DataType> > newField(new Field<Dimension, DataType>(**itr));
      mFieldCache.push_back(newField);
      *itr = mFieldCache.back().get();
      *baseItr = mFieldCache.back().get();
    }

//     for (int i = 0; i < this->size(); ++i) {
//       mFieldCache.push_back(*((*this)[i]));
//       (*this)[i] = &(mFieldCache.back());
//     }
  }
}

//------------------------------------------------------------------------------
// Test if a given field is part of this field list.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
haveField(const Field<Dimension, DataType>& field) const {
  const_iterator fieldListItr = std::find(this->begin(), this->end(), &field);
  return fieldListItr != this->end();
}

//------------------------------------------------------------------------------
// Test if there is a Field associated with the given NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
haveNodeList(const NodeList<Dimension>& nodeList) const {
  return mNodeListIndexMap.find(&nodeList) != mNodeListIndexMap.end();
}

//------------------------------------------------------------------------------
// Set the fields of this field list equal to those of another.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
assignFields(const FieldList<Dimension, DataType>& fieldList) {
  const_iterator otherFieldListItr = fieldList.begin();
  for (iterator fieldListItr = this->begin(); fieldListItr < this->end();
       ++fieldListItr, ++otherFieldListItr) {
    CHECK(otherFieldListItr < fieldList.end());
    CHECK((*fieldListItr)->nodeListPtr() == (*otherFieldListItr)->nodeListPtr());
    **fieldListItr = **otherFieldListItr;
  }
}

//------------------------------------------------------------------------------
// Append the given field to the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::appendField(const Field<Dimension, DataType>& field) {
  if (haveField(field)) {
    std::cerr << "FieldList::appendField Warning: attempt to append field " << &field
              << " to FieldList " << this
              << " which already has it." << std::endl;
    return;
  }

  // Determine the order this Field should be in.
  const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
  iterator orderItr = nlr.findInsertionPoint(&field,
                                             begin(),
                                             end());
  const int delta = std::distance(begin(), orderItr);

  // Insert the field.
  switch(storageType()) {
  case FieldStorageType::ReferenceFields:
    mFieldPtrs.insert(orderItr, const_cast<Field<Dimension, DataType>*>(&field));
    mFieldBasePtrs.insert(mFieldBasePtrs.begin() + delta, const_cast<FieldBase<Dimension>*>(dynamic_cast<const FieldBase<Dimension>*>(&field)));
    break;

  case FieldStorageType::CopyFields:
    std::shared_ptr<Field<Dimension, DataType> > newField(new Field<Dimension, DataType>(field));
    mFieldCache.push_back(newField);
    mFieldPtrs.insert(orderItr, newField.get());
    mFieldBasePtrs.insert(mFieldBasePtrs.begin() + delta, newField.get());
  }

//   registerWithField(*mFieldPtrs.back());

  // We also update the set of NodeListPtrs in proper order.
  mNodeListPtrs.insert(mNodeListPtrs.begin() + delta,
                       const_cast<NodeList<Dimension>*>(field.nodeListPtr()));
  CHECK(mNodeListIndexMap.find(field.nodeListPtr()) == mNodeListIndexMap.end());
  buildNodeListIndexMap();

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  ENSURE(mFieldPtrs[mNodeListIndexMap[field.nodeListPtr()]] == *(fieldForNodeList(*field.nodeListPtr())));
  ENSURE(this->size() == mNodeListPtrs.size());
  END_CONTRACT_SCOPE

}

//------------------------------------------------------------------------------
// Delete the given field from the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::deleteField(const Field<Dimension, DataType>& field) {
  if (!haveField(field)) {
    std::cerr << "FieldList::deleteField Warning: attempt to delete field " << &field
              << " from FieldList " << this
              << " which does not recognize it." << std::endl;
    return;
  }

  const iterator fieldPtrItr = std::find(this->begin(), this->end(), &field);
  CHECK(fieldPtrItr != this->end());
  const size_t delta = std::distance(this->begin(), fieldPtrItr);
  typename FieldCacheType::iterator fieldItr = mFieldCache.begin();
  switch(storageType()) {
  case FieldStorageType::CopyFields:
    while (fieldItr != mFieldCache.end() && fieldItr->get() != &field) {
      ++fieldItr;
    }
    CHECK(fieldItr != mFieldCache.end());
    mFieldCache.erase(fieldItr);

  case FieldStorageType::ReferenceFields:
    mFieldPtrs.erase(fieldPtrItr);
    mFieldBasePtrs.erase(mFieldBasePtrs.begin() + delta);
    break;
  }

  // Remove the NodeList pointer.
  typename std::vector<NodeList<Dimension>*>::iterator nodeListItr = std::find(mNodeListPtrs.begin(), mNodeListPtrs.end(), field.nodeListPtr());
  CHECK(nodeListItr != mNodeListPtrs.end());
  mNodeListPtrs.erase(nodeListItr);
  buildNodeListIndexMap();
}

//------------------------------------------------------------------------------
// Construct a new field and add it to the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
appendNewField(const typename Field<Dimension, DataType>::FieldName name,
               const NodeList<Dimension>& nodeList,
               const DataType value) {
  VERIFY(mStorageType == FieldStorageType::CopyFields);

  // Create the field in our cache.
  mFieldCache.push_back(std::shared_ptr<Field<Dimension, DataType> >(new Field<Dimension, DataType>(name, nodeList, value)));
  Field<Dimension, DataType>* fieldPtr = mFieldCache.back().get();

  // Determine the order this Field should be in.
  const NodeListRegistrar<Dimension>& nlr = NodeListRegistrar<Dimension>::instance();
  iterator orderItr = nlr.findInsertionPoint(fieldPtr,
                                             begin(),
                                             end());
  const int delta = std::distance(begin(), orderItr);

  // Insert the field.
  mFieldPtrs.insert(orderItr, fieldPtr);
  mFieldBasePtrs.insert(mFieldBasePtrs.begin() + delta, fieldPtr);

  // We also update the set of NodeListPtrs in proper order.
  mNodeListPtrs.insert(mNodeListPtrs.begin() + delta, const_cast<NodeList<Dimension>*>(&nodeList));
  CHECK(mNodeListIndexMap.find(fieldPtr->nodeListPtr()) == mNodeListIndexMap.end());
  buildNodeListIndexMap();

  // Post-conditions.
  BEGIN_CONTRACT_SCOPE
  ENSURE(mFieldPtrs[mNodeListIndexMap[fieldPtr->nodeListPtr()]] == *(fieldForNodeList(*(fieldPtr->nodeListPtr()))));
  ENSURE(this->size() == mNodeListPtrs.size());
  ENSURE(mFieldBasePtrs.size() == mFieldPtrs.size());
  END_CONTRACT_SCOPE

}

//------------------------------------------------------------------------------
// Standard iterators.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::iterator
FieldList<Dimension, DataType>::
begin() {
  return mFieldPtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::iterator
FieldList<Dimension, DataType>::
end() {
  return mFieldPtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::reverse_iterator
FieldList<Dimension, DataType>::
rbegin() {
  return mFieldPtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::reverse_iterator
FieldList<Dimension, DataType>::
rend() {
  return mFieldPtrs.rend();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_iterator
FieldList<Dimension, DataType>::
begin() const {
  return mFieldPtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_iterator
FieldList<Dimension, DataType>::
end() const {
  return mFieldPtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_reverse_iterator
FieldList<Dimension, DataType>::
rbegin() const {
  return mFieldPtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_reverse_iterator
FieldList<Dimension, DataType>::
rend() const {
  return mFieldPtrs.rend();
}

//------------------------------------------------------------------------------
// Standard iterators for FieldListBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::iterator
FieldList<Dimension, DataType>::
begin_base() {
  return mFieldBasePtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::iterator
FieldList<Dimension, DataType>::
end_base() {
  return mFieldBasePtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::reverse_iterator
FieldList<Dimension, DataType>::
rbegin_base() {
  return mFieldBasePtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::reverse_iterator
FieldList<Dimension, DataType>::
rend_base() {
  return mFieldBasePtrs.rend();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_iterator
FieldList<Dimension, DataType>::
begin_base() const {
  return mFieldBasePtrs.begin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_iterator
FieldList<Dimension, DataType>::
end_base() const {
  return mFieldBasePtrs.end();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_reverse_iterator
FieldList<Dimension, DataType>::
rbegin_base() const {
  return mFieldBasePtrs.rbegin();
}

template<typename Dimension, typename DataType>
inline
typename FieldListBase<Dimension>::const_reverse_iterator
FieldList<Dimension, DataType>::
rend_base() const {
  return mFieldBasePtrs.rend();
}

//------------------------------------------------------------------------------
// Index operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
operator[](const unsigned int index) const {
  REQUIRE(index < mFieldPtrs.size());
  return mFieldPtrs[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
operator[](const unsigned int index) {
  REQUIRE(index < mFieldPtrs.size());
  return mFieldPtrs[index];
}

//------------------------------------------------------------------------------
// at version, for consistency with the STL.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
at(const unsigned int index) const {
  return (*this)[index];
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::ElementType
FieldList<Dimension, DataType>::
at(const unsigned int index) {
  return (*this)[index];
}

//------------------------------------------------------------------------------
// Return an iterator pointing to the Field member associated with the given
// NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::iterator
FieldList<Dimension, DataType>::
fieldForNodeList(const NodeList<Dimension>& nodeList) {
  if (haveNodeList(nodeList)) {
    return begin() + mNodeListIndexMap.find(&nodeList)->second;
  } else {
    return end();
  }
}

template<typename Dimension, typename DataType>
inline
typename FieldList<Dimension, DataType>::const_iterator
FieldList<Dimension, DataType>::
fieldForNodeList(const NodeList<Dimension>& nodeList) const {
  if (haveNodeList(nodeList)) {
    return begin() + mNodeListIndexMap.find(&nodeList)->second;
  } else {
    return end();
  }
}

//------------------------------------------------------------------------------
// Element access by NodeIteratorBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldList<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) {
  return this->operator()(itr.fieldID(), itr.nodeID());
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldList<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) const {
  return this->operator()(itr.fieldID(), itr.nodeID());
}

//------------------------------------------------------------------------------
// Provide a more primitive access to Field elements, based on the index of the Field
// and the node index within that field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
FieldList<Dimension, DataType>::
operator()(const unsigned int fieldIndex,
           const unsigned int nodeIndex) {
  REQUIRE(fieldIndex < mFieldPtrs.size());
  REQUIRE(nodeIndex < mFieldPtrs[fieldIndex]->size());
  return mFieldPtrs[fieldIndex]->operator()(nodeIndex);
}

template<typename Dimension, typename DataType>
inline
const DataType&
FieldList<Dimension, DataType>::
operator()(const unsigned int fieldIndex,
           const unsigned int nodeIndex) const {
  REQUIRE(fieldIndex < mFieldPtrs.size());
  REQUIRE(nodeIndex < mFieldPtrs[fieldIndex]->size());
  return mFieldPtrs[fieldIndex]->operator()(nodeIndex);
}

//------------------------------------------------------------------------------
// Return the interpolated value of the FieldList for the given position.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
operator()(const typename Dimension::Vector& position,
           const TableKernel<Dimension>& W) const {

  DataType result(0.0);

  // Set the neighbor information for all NodeLists.
  std::vector<std::vector<int>> masterLists, coarseNeighbors, refineNeighbors;
  setMasterNodeLists(position, masterLists, coarseNeighbors);
  setRefineNodeLists(position, coarseNeighbors, refineNeighbors);

  // Loop over the neighbors.
  for (RefineNodeIterator<Dimension> neighborItr = refineNodeBegin(refineNeighbors);
       neighborItr < refineNodeEnd();
       ++neighborItr) {

    // Ignore any ghost nodes.
    if (neighborItr.internalNode()) {

      // This is horrible, but I need the node weights to do this job.  This routine
      // needs to be rewritten!
      const Scalar& wj = (neighborItr.fluidNodeListPtr()->weight())(neighborItr);

      // Add this nodes contribution.
      const Vector rij = position - neighborItr.fluidNodeListPtr()->positions()(neighborItr);
      const SymTensor& Hj = neighborItr.fluidNodeListPtr()->Hfield()(neighborItr);
      double Wj = W(Hj*rij, Hj);
      result += (*this)(neighborItr)*wj*Wj;
    }
  }

#ifdef USE_MPI
  // In parallel, we need to sum up the result across all processors.
  {
    int procID, numProcs;
    MPI_Comm_rank(Communicator::communicator(), &procID);
    MPI_Comm_size(Communicator::communicator(), &numProcs);
    VERIFY(DataTypeTraits<DataType>::fixedSize());
    const size_t sizeOfElement = DataTypeTraits<DataType>::numElements()*sizeof(typename DataTypeTraits<DataType>::ElementType);
    std::vector<char> sendBuffer;
    std::vector<char> recvBuffer(sizeOfElement);
    packElement(result, sendBuffer);
    CHECK(sendBuffer.size() == sizeOfElement);
    for (int sendProc = 0; sendProc != numProcs; ++sendProc) {
      recvBuffer = sendBuffer;
      MPI_Bcast(&(*recvBuffer.begin()), sizeOfElement, MPI_CHAR, sendProc, Communicator::communicator());
      if (procID != sendProc) {
        DataType otherValue;
        std::vector<char>::const_iterator itr = recvBuffer.begin();
        unpackElement(otherValue, itr, recvBuffer.end());
        CHECK(itr == recvBuffer.end());
        result += otherValue;
      }
    }
  }
#endif

  return result;
}

//------------------------------------------------------------------------------
// Node ID iterators for all node IDs.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
AllNodeIterator<Dimension>
FieldList<Dimension, DataType>::nodeBegin() const {
  typename std::vector<NodeList<Dimension>*>::const_iterator 
    nodeListItr = mNodeListPtrs.begin();
  while (nodeListItr < mNodeListPtrs.end() &&
         (*nodeListItr)->numNodes() == 0) {
    ++nodeListItr;
  }
  return AllNodeIterator<Dimension>(nodeListItr,
                                    mNodeListPtrs.begin(),
                                    mNodeListPtrs.end());
}

template<typename Dimension, typename DataType>
inline
AllNodeIterator<Dimension>
FieldList<Dimension, DataType>::nodeEnd() const {
  return AllNodeIterator<Dimension>(mNodeListPtrs.end(),
                                    mNodeListPtrs.begin(),
                                    mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Node ID iterators for internal nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
InternalNodeIterator<Dimension>
FieldList<Dimension, DataType>::internalNodeBegin() const {
  typename std::vector<NodeList<Dimension>*>::const_iterator
    nodeListItr = mNodeListPtrs.begin();
  while (nodeListItr < mNodeListPtrs.end() &&
         (*nodeListItr)->numInternalNodes() == 0) {
    ++nodeListItr;
  }
  return InternalNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end());
}

template<typename Dimension, typename DataType>
inline
InternalNodeIterator<Dimension>
FieldList<Dimension, DataType>::internalNodeEnd() const {
  return InternalNodeIterator<Dimension>(mNodeListPtrs.end(),
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Node ID iterators for ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
GhostNodeIterator<Dimension>
FieldList<Dimension, DataType>::ghostNodeBegin() const {
  typename std::vector<NodeList<Dimension>*>::const_iterator
    nodeListItr = mNodeListPtrs.begin();
  while (nodeListItr < mNodeListPtrs.end() &&
         (*nodeListItr)->numGhostNodes() == 0) {
    ++nodeListItr;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    CHECK((*nodeListItr)->firstGhostNode() < (*nodeListItr)->numNodes());
    return GhostNodeIterator<Dimension>(nodeListItr,
                                        mNodeListPtrs.begin(),
                                        mNodeListPtrs.end(),
                                        (*nodeListItr)->firstGhostNode());
  } else {
    return ghostNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
GhostNodeIterator<Dimension>
FieldList<Dimension, DataType>::ghostNodeEnd() const {
  return GhostNodeIterator<Dimension>(mNodeListPtrs.end(),
                                      mNodeListPtrs.begin(),
                                      mNodeListPtrs.end());
}

//------------------------------------------------------------------------------
// Node ID iterators for master neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
MasterNodeIterator<Dimension>
FieldList<Dimension, DataType>::masterNodeBegin(const std::vector<std::vector<int>>& masterLists) const {
  auto nodeListItr = mNodeListPtrs.begin();
  unsigned iNodeList = 0;
  while (nodeListItr < mNodeListPtrs.end() && masterLists[iNodeList].empty()) {
    ++nodeListItr;
    ++iNodeList;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    return MasterNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end(),
                                         masterLists[iNodeList].begin(),
                                         masterLists);
  } else {
    return this->masterNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
MasterNodeIterator<Dimension>
FieldList<Dimension, DataType>::masterNodeEnd() const {
  return MasterNodeIterator<Dimension>(mNodeListPtrs.end(),
                                       mNodeListPtrs.begin(),
                                       mNodeListPtrs.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node ID iterators for coarse neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
CoarseNodeIterator<Dimension>
FieldList<Dimension, DataType>::coarseNodeBegin(const std::vector<std::vector<int>>& coarseNeighbors) const {
  auto nodeListItr = mNodeListPtrs.begin();
  unsigned iNodeList = 0;
  while (nodeListItr < mNodeListPtrs.end() && coarseNeighbors[iNodeList].empty()) {
    ++nodeListItr;
    ++iNodeList;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    return CoarseNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end(),
                                         coarseNeighbors[iNodeList].begin(),
                                         coarseNeighbors);
  } else {
    return this->coarseNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
CoarseNodeIterator<Dimension>
FieldList<Dimension, DataType>::coarseNodeEnd() const {
  return CoarseNodeIterator<Dimension>(mNodeListPtrs.end(),
                                       mNodeListPtrs.begin(),
                                       mNodeListPtrs.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Node ID iterators for refine neighbor nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
RefineNodeIterator<Dimension>
FieldList<Dimension, DataType>::refineNodeBegin(const std::vector<std::vector<int>>& refineNeighbors) const {
  auto nodeListItr = mNodeListPtrs.begin();
  unsigned iNodeList = 0;
  while (nodeListItr < mNodeListPtrs.end() && refineNeighbors[iNodeList].empty()) {
    ++nodeListItr;
    ++iNodeList;
  }
  if (nodeListItr < mNodeListPtrs.end()) {
    return RefineNodeIterator<Dimension>(nodeListItr,
                                         mNodeListPtrs.begin(),
                                         mNodeListPtrs.end(),
                                         refineNeighbors[iNodeList].begin(),
                                         refineNeighbors);
  } else {
    return this->refineNodeEnd();
  }
}

template<typename Dimension, typename DataType>
inline
RefineNodeIterator<Dimension>
FieldList<Dimension, DataType>::refineNodeEnd() const {
  return RefineNodeIterator<Dimension>(mNodeListPtrs.end(),
                                       mNodeListPtrs.begin(),
                                       mNodeListPtrs.end(),
                                       std::vector<std::vector<int>>());
}

//------------------------------------------------------------------------------
// Set the master node lists of all the NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setMasterNodeLists(const typename Dimension::Vector& r,
                   const typename Dimension::SymTensor& H,
                   std::vector<std::vector<int>>& masterLists,
                   std::vector<std::vector<int>>& coarseNeighbors) const {
  auto etaMax = 0.0;
  for (auto nodeListItr = mNodeListPtrs.begin();
       nodeListItr != mNodeListPtrs.end();
       ++nodeListItr) etaMax = std::max(etaMax, (**nodeListItr).neighbor().kernelExtent());
  Neighbor<Dimension>::setMasterNeighborGroup(r, H,
                                                             mNodeListPtrs.begin(),
                                                             mNodeListPtrs.end(),
                                                             etaMax,
                                                             masterLists,
                                                             coarseNeighbors);
}

template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setMasterNodeLists(const typename Dimension::Vector& r,
                   std::vector<std::vector<int>>& masterLists,
                   std::vector<std::vector<int>>& coarseNeighbors) const {
  this->setMasterNodeLists(r, 1e-30*SymTensor::one, masterLists, coarseNeighbors);
}

//------------------------------------------------------------------------------
// Set the refine node lists of all the NodeLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setRefineNodeLists(const typename Dimension::Vector& r,
                   const typename Dimension::SymTensor& H,
                   const std::vector<std::vector<int>>& coarseNeighbors,
                   std::vector<std::vector<int>>& refineNeighbors) const {
  const auto numNodeLists = mNodeListPtrs.size();
  REQUIRE(coarseNeighbors.size() == numNodeLists);
  refineNeighbors = std::vector<std::vector<int>>(numNodeLists);
  auto iNodeList = 0;
  for (auto nodeListItr = mNodeListPtrs.begin();
       nodeListItr < mNodeListPtrs.end();
       ++nodeListItr, ++iNodeList) {
    (*nodeListItr)->neighbor().setRefineNeighborList(r, H, 
                                                     coarseNeighbors[iNodeList],
                                                     refineNeighbors[iNodeList]);
  }
}

template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
setRefineNodeLists(const typename Dimension::Vector& r,
                   const std::vector<std::vector<int>>& coarseNeighbors,
                   std::vector<std::vector<int>>& refineNeighbors) const {
  this->setRefineNodeLists(r, 1e-30*SymTensor::one, coarseNeighbors, refineNeighbors);
}

//------------------------------------------------------------------------------
// Zero out the FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::Zero() {
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->Zero();
  }
}

//------------------------------------------------------------------------------
// Apply a minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyMin(const DataType& dataMin) {
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyMin(dataMin);
  }
}

//------------------------------------------------------------------------------
// Apply a maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyMax(const DataType& dataMax) {
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyMax(dataMax);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) minimum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyScalarMin(const double dataMin) {
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyScalarMin(dataMin);
  }
}

//------------------------------------------------------------------------------
// Apply a (double) maximum data value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::applyScalarMax(const double dataMax) {
  for (iterator fieldItr = begin(); fieldItr < end(); ++fieldItr) {
    (*fieldItr)->applyScalarMax(dataMax);
  }
}

//------------------------------------------------------------------------------
// Add two FieldLists.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator+(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a FieldList from another.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator-(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Add two FieldLists in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator+=(const FieldList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  for (int i = 0; i < numFields(); ++i) {
    CHECK((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr());
    *((*this)[i]) += *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a FieldList from another in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator-=(const FieldList<Dimension, DataType>& rhs) {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  for (int i = 0; i < numFields(); ++i) {
    CHECK((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr());
    *((*this)[i]) -= *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Add a single value to a FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator+(const DataType& rhs) const {
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a single value from a Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator-(const DataType& rhs) const {
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Add a single value to the Field in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator+=(const DataType& rhs) {
  for (int i = 0; i < numFields(); ++i) {
    *((*this)[i]) += rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Subtract a single value from the FieldList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator-=(const DataType& rhs) {
  for (int i = 0; i < numFields(); ++i) {
    *((*this)[i]) -= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldList by a Scalar FieldList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator*=(const FieldList<Dimension, typename Dimension::Scalar>& rhs) {
  REQUIRE(this->numFields() == rhs.numFields());
  for (int i = 0; i < numFields(); ++i) {
    CHECK((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr());
    *((*this)[i]) *= *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Multiply this FieldList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator*=(const Scalar& rhs) {
  for (int i = 0; i < numFields(); ++i) {
    *((*this)[i]) *= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::
operator/(const FieldList<Dimension, typename Dimension::Scalar>& rhs) const {
  REQUIRE(this->numFields() == rhs.numFields());
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>
FieldList<Dimension, DataType>::operator/(const Scalar& rhs) const {
  REQUIRE(rhs != 0.0);
  FieldList<Dimension, DataType> result(*this);
  result.copyFields();
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar FieldList in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::
operator/=(const FieldList<Dimension, typename Dimension::Scalar>& rhs) {
  REQUIRE(this->numFields() == rhs.numFields());
  for (int i = 0; i < numFields(); ++i) {
    CHECK((*this)[i]->nodeListPtr() == rhs[i]->nodeListPtr());
    *((*this)[i]) /= *(rhs[i]);
  }
  return *this;
}

//------------------------------------------------------------------------------
// Divide this FieldList by a Scalar in place.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldList<Dimension, DataType>&
FieldList<Dimension, DataType>::operator/=(const typename Dimension::Scalar& rhs) {
  REQUIRE(rhs != 0.0);
  for (int i = 0; i < numFields(); ++i) {
    *((*this)[i]) /= rhs;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Sum the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
sumElements() const {
  DataType result = DataTypeTraits<DataType>::zero();
  for (const_iterator itr = begin(); itr != end(); ++itr) {
    result += (*itr)->sumElements();
  }
  return result;
}

//------------------------------------------------------------------------------
// Find the minimum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
min() const {
  DataType result = std::numeric_limits<DataType>::max();
  for (const_iterator itr = begin(); itr != end(); ++itr) {
    result = std::min(result, (*itr)->min());
  }
  return result;
}

//------------------------------------------------------------------------------
// Find the maximum.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
max() const {
  DataType result;
  result = -std::numeric_limits<DataType>::max() < std::numeric_limits<DataType>::min() ? -std::numeric_limits<DataType>::max() : std::numeric_limits<DataType>::min();
  for (const_iterator itr = begin(); itr != end(); ++itr) {
    result = std::max(result, (*itr)->max());
  }
  return result;
}

//------------------------------------------------------------------------------
// Sum the field elements.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
localSumElements() const {
  DataType result = DataTypeTraits<DataType>::zero();
  for (const_iterator itr = begin(); itr != end(); ++itr) {
    result += (*itr)->localSumElements();
  }
  return result;
}

//------------------------------------------------------------------------------
// Find the minimum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
localMin() const {
  DataType result = std::numeric_limits<DataType>::max();
  for (const_iterator itr = begin(); itr != end(); ++itr) {
    result = std::min(result, (*itr)->localMin());
  }
  return result;
}

//------------------------------------------------------------------------------
// Find the maximum.
// LOCAL to processor!
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType
FieldList<Dimension, DataType>::
localMax() const {
  DataType result;
  result = -std::numeric_limits<DataType>::max() < std::numeric_limits<DataType>::min() ? -std::numeric_limits<DataType>::max() : std::numeric_limits<DataType>::min();
  for (const_iterator itr = begin(); itr != end(); ++itr) {
    result = std::max(result, (*itr)->localMax());
  }
  return result;
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator==(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) == *(rhs[i]));
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator!=(const FieldList<Dimension, DataType>& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) > *(rhs[i]));
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<(const FieldList<Dimension, DataType>& rhs) const {

  // Pre-conditions.
  BEGIN_CONTRACT_SCOPE
  REQUIRE(numFields() == rhs.numFields());
  for (int i = 0; i != numFields(); ++i) 
    REQUIRE(mFieldPtrs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
  END_CONTRACT_SCOPE

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) < *(rhs[i]));
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>=(const FieldList<Dimension, DataType>& rhs) const {
  return operator>(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<=(const FieldList<Dimension, DataType>& rhs) const {
  return operator<(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator==
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator==(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) == rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator!=(const DataType& rhs) const {
  return !(operator==(rhs));
}

//------------------------------------------------------------------------------
// operator>
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) > rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator<
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<(const DataType& rhs) const {

  bool result = true;
  int i = 0;
  while (result && i != numFields()) {
    result = result && (*(mFieldPtrs[i]) < rhs);
    ++i;
  }
  return result;
}

//------------------------------------------------------------------------------
// operator>=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator>=(const DataType& rhs) const {
  return operator>(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// operator<=
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
FieldList<Dimension, DataType>::
operator<=(const DataType& rhs) const {
  return operator<(rhs) || operator==(rhs);
}

//------------------------------------------------------------------------------
// Return the number of Fields stored in this Field List.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
unsigned 
FieldList<Dimension, DataType>::numFields() const {
  return mFieldPtrs.size();
}

template<typename Dimension, typename DataType>
inline
unsigned 
FieldList<Dimension, DataType>::size() const {
  return numFields();
}

template<typename Dimension, typename DataType>
inline
unsigned 
FieldList<Dimension, DataType>::numNodes() const {
  unsigned numberOfNodes = 0;
  for (const_iterator iter = begin(); iter != end(); ++iter) {
    numberOfNodes += (*iter)->nodeList().numNodes();
  }
  return numberOfNodes;
}

template<typename Dimension, typename DataType>
inline
unsigned 
FieldList<Dimension, DataType>::numInternalNodes() const {
  unsigned numberOfNodes = 0;
  for (const_iterator iter = begin(); iter != end(); ++iter) {
    numberOfNodes += (*iter)->nodeList().numInternalNodes();
  }
  return numberOfNodes;
}

template<typename Dimension, typename DataType>
inline
unsigned 
FieldList<Dimension, DataType>::numGhostNodes() const {
  unsigned numberOfNodes = 0;
  for (const_iterator iter = begin(); iter != end(); ++iter) {
    numberOfNodes += (*iter)->nodeList().numGhostNodes();
  }
  return numberOfNodes;
}

//------------------------------------------------------------------------------
// Return the set of NodeList pointers.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
const std::vector<NodeList<Dimension>*>&
FieldList<Dimension, DataType>::
nodeListPtrs() const {
  return mNodeListPtrs;
}

//------------------------------------------------------------------------------
// Internal method to build the NodeListIndexMap from scratch.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
FieldList<Dimension, DataType>::
buildNodeListIndexMap() {
  mNodeListIndexMap = HashMapType();
  int i = 0;
  for (iterator itr = begin();
       itr != end();
       ++itr, ++i) mNodeListIndexMap[(*itr)->nodeListPtr()] = i;
  ENSURE(mNodeListIndexMap.size() == numFields());
}

//****************************** Global Functions ******************************

//------------------------------------------------------------------------------
// Multiply a FieldList by another FieldList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
inline
FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const FieldList<Dimension, DataType>& lhs,
          const FieldList<Dimension, OtherDataType>& rhs) {
  REQUIRE(lhs.numFields() == rhs.numFields());
  FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
  result.copyFields();
  for (int i = 0; i < lhs.numFields(); ++i) {
    CHECK(lhs[i]->nodeListPtr() == rhs[i]->nodeListPtr());
    result.appendField((*(lhs[i])) * (*(rhs[i])));
  }
  return result;
}

//------------------------------------------------------------------------------
// Multiply a FieldList by a single value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
inline
FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const FieldList<Dimension, DataType>& lhs,
          const OtherDataType& rhs) {
  FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
  result.copyFields();
  for (int i = 0; i < lhs.numFields(); ++i) {
    result.appendField((*(lhs[i])) * rhs);
  }
  return result;
}

template<typename Dimension, typename DataType, typename OtherDataType>
inline
FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const DataType& lhs,
          const FieldList<Dimension, OtherDataType>& rhs) {
  FieldList<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType> result;
  result.copyFields();
  for (int i = 0; i < rhs.numFields(); ++i) {
    result.appendField(lhs * (*(rhs[i])));
  }
  return result;
}

// //------------------------------------------------------------------------------
// // Absolute value.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// abs(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(abs(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse cosine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// acos(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(acos(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse sine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// asin(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(asin(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse tangent.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// atan(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(atan(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Inverse tangent2.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// atan2(const FieldList<Dimension, typename Dimension::Scalar>& fieldList1,
//       const FieldList<Dimension, typename Dimension::Scalar>& fieldList2) {
//   REQUIRE(fieldList1.numFields() == fieldList2.numFields());
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList1.numFields(); ++i) {
//     result.appendField(atan2(*(fieldList1[i]), *(fieldList2[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Ceiling.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// ceil(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(ceil(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Cosine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// cos(FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(cos(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Hyperbolic cosine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// cosh(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(cosh(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Exponential.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// exp(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(exp(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Absolute value.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// fabs(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(fabs(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Floor.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// floor(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(floor(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Natural logarithm.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// log(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(log(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Log base 10.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// log10(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(log10(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // powN.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow2(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow2(*(fieldList[i])));
//   }
//   return result;
// }


// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow3(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow3(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow4(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow4(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow5(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow5(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow6(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow6(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow7(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow7(*(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// pow8(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(pow8(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Sine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sin(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sin(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Hyperbolic sine.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sinh(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sinh(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Square.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sqr(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sqr(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Square root.
// //------------------------------------------------------------------------------
// template<typename Dimension>
// inline
// FieldList<Dimension, typename Dimension::Scalar>
// sqrt(const FieldList<Dimension, typename Dimension::Scalar>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(sqrt(*(fieldList[i])));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Minimum.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// min(const FieldList<Dimension, DataType>& fieldList1,
//     const FieldList<Dimension, DataType>& fieldList2) {
//   REQUIRE(fieldList1.numFields() == fieldList2.numFields());
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList1.numFields(); ++i) {
//     result.appendField(min(*(fieldList1[i]), *(fieldList2[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// min(const DataType& value,
//     const FieldList<Dimension, DataType>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(min(value, *(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// min(const FieldList<Dimension, DataType>& fieldList,
//     const DataType& value) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(min(*(fieldList[i]), value));
//   }
//   return result;
// }

// //------------------------------------------------------------------------------
// // Maximum.
// //------------------------------------------------------------------------------
// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// max(const FieldList<Dimension, DataType>& fieldList1,
//     const FieldList<Dimension, DataType>& fieldList2) {
//   REQUIRE(fieldList1.numFields() == fieldList2.numFields());
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList1.numFields(); ++i) {
//     result.appendField(max(*(fieldList1[i]), *(fieldList2[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// max(const DataType& value,
//     const FieldList<Dimension, DataType>& fieldList) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(max(value, *(fieldList[i])));
//   }
//   return result;
// }

// template<typename Dimension, typename DataType>
// inline
// FieldList<Dimension, DataType>
// max(const FieldList<Dimension, DataType>& fieldList,
//     const DataType& value) {
//   FieldList<Dimension, typename Dimension::Scalar> result;
//   result.copyFields();
//   for (int i = 0; i < fieldList.numFields(); ++i) {
//     result.appendField(max(*(fieldList[i]), value));
//   }
//   return result;
// }

}
