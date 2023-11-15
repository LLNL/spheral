#include "Field/FieldBase.hh"
#include "Geometry/Dimension.hh"
#include "NodeList/NodeList.hh"
#include "Field/NodeIterators.hh"
#include "Utilities/packElement.hh"
#include "Utilities/removeElements.hh"
#include "Utilities/safeInv.hh"
#include "Utilities/allReduce.hh"
#include "Distributed/Communicator.hh"

#include <cmath>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <sstream>
#include <limits>

#ifdef USE_MPI
extern "C" {
#include <mpi.h>
}
#endif

// Inlined methods.
namespace Spheral {

// Construct with name.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name):
  FieldBase<Dimension>(name),
  mValid(false) {}

////------------------------------------------------------------------------------
//// Construct with name and field values.
////------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const Field<Dimension, DataType>& field):
  FieldBase<Dimension>(name, *field.nodeListPtr()),
  mValid(field.mValid) {
    FieldViewType::mDataArray = ContainerType(deepCopy(field.mDataArray));
  }

//------------------------------------------------------------------------------
// Construct with the given name and NodeList.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const NodeList<Dimension>& nodeList):
  FieldBase<Dimension>(name, nodeList),
  mValid(true) {
  FieldViewType::mDataArray = ContainerType((size_t) nodeList.numNodes(), DataType());
  REQUIRE(numElements() == nodeList.numNodes());
}

template<>
inline
Field<Dim<1>, Dim<1>::Scalar>::
Field(FieldBase<Dim<1> >::FieldName name,
      const NodeList<Dim<1> >& nodeList):
  FieldBase<Dim<1> >(name, nodeList),
  mValid(true) {
  FieldViewType::mDataArray = ContainerType((size_t) nodeList.numNodes(), 0.0);
  REQUIRE(numElements() == nodeList.numNodes());
}

template<>
inline
Field<Dim<2>, Dim<2>::Scalar>::
Field(FieldBase<Dim<2> >::FieldName name,
      const NodeList<Dim<2> >& nodeList):
  FieldBase<Dim<2> >(name, nodeList),
  mValid(true) {
  FieldViewType::mDataArray = ContainerType((size_t) nodeList.numNodes(), 0.0);
  REQUIRE(numElements() == nodeList.numNodes());
}

template<>
inline
Field<Dim<3>, Dim<3>::Scalar>::
Field(FieldBase<Dim<3> >::FieldName name,
      const NodeList<Dim<3> >& nodeList):
  FieldBase<Dim<3> >(name, nodeList),
  mValid(true) {
  FieldViewType::mDataArray = ContainerType((size_t) nodeList.numNodes(), 0.0);
  REQUIRE(numElements() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Construct with given name, NodeList, and value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name,
      const NodeList<Dimension>& nodeList,
      DataType value):
  FieldBase<Dimension>(name, nodeList),
  mValid(true) {
  FieldViewType::mDataArray = ContainerType((size_t) nodeList.numNodes(), value);
  REQUIRE(numElements() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Construct for a given name, NodeList, and vector of values.
//------------------------------------------------------------------------------

template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::
Field(typename FieldBase<Dimension>::FieldName name, 
      const NodeList<Dimension>& nodeList,
      const ContainerType& array):
  FieldBase<Dimension>(name, nodeList),
  mValid(true) {
  REQUIRE(numElements() == nodeList.numNodes());
  REQUIRE(numElements() == array.size());
  FieldViewType::mDataArray = ContainerType(deepCopy(array));
}

//------------------------------------------------------------------------------
// Construct by copying the values of another Field, but using a different
// node list.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::Field(const NodeList<Dimension>& nodeList,
                                  const Field<Dimension, DataType>& field):
  FieldBase<Dimension>(field.name(), nodeList),
  mValid(true) {
  FieldViewType::mDataArray = ContainerType(deepCopy(field.mDataArray));
  ENSURE(numElements() == nodeList.numNodes());
}

//------------------------------------------------------------------------------
// Copy Constructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::Field(const Field& field):
  FieldBase<Dimension>(field),
  FieldView<Dimension, DataType>(field),
  mValid(field.valid()) {
  FieldViewType::mDataArray = ContainerType(deepCopy(field.mDataArray));
}

//------------------------------------------------------------------------------
// The virtual clone method, allowing us to duplicate fields with just 
// FieldBase*.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::shared_ptr<FieldBase<Dimension> >
Field<Dimension, DataType>::clone() const {
  return std::shared_ptr<FieldBase<Dimension>>(new Field<Dimension, DataType>(*this));
}

//------------------------------------------------------------------------------
// Destructor.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>::~Field() {
  FieldViewType::mDataArray.free();
}

//------------------------------------------------------------------------------
// Assignment operator with FieldBase.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
FieldBase<Dimension>&
Field<Dimension, DataType>::operator=(const FieldBase<Dimension>& rhs) {
  if (this != &rhs) {
    //TODO: Figure out a way to handle this on the device later...
#if !defined(SPHERAL_GPU_ACTIVE)
    try {
      const Field<Dimension, DataType>* rhsPtr = dynamic_cast<const Field<Dimension, DataType>*>(&rhs);
      CHECK2(rhsPtr != 0, "Passed incorrect Field to operator=!");
      FieldBase<Dimension>::operator=(rhs);
      FieldViewType::mDataArray = deepCopy(rhsPtr->mDataArray);
      mValid = rhsPtr->mValid;
    } catch (const std::bad_cast &) {
      VERIFY2(false, "Attempt to assign a field to an incompatible field type.");
    }
#endif // SPHERAL_GPU_ACTIVE
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assignment operator to another Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const Field<Dimension, DataType>& rhs) {
  REQUIRE(rhs.valid());
  if (this != &rhs) {
    FieldBase<Dimension>::operator=(rhs);
    FieldViewType::mDataArray = deepCopy(rhs.mDataArray);
    mValid = rhs.mValid;
  }
  return *this;
}

//------------------------------------------------------------------------------
// Assigment operator with a vector.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const ContainerType& rhs) {
  REQUIRE(mValid);
  REQUIRE(this->nodeList().numNodes() == rhs.size());
  FieldViewType::mDataArray = deepCopy(rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Assignment operator with a constant value of DataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>&
Field<Dimension, DataType>::operator=(const DataType& rhs) {
  REQUIRE(mValid);
  std::fill(FieldViewType::begin(), FieldViewType::end(), rhs);
  return *this;
}

//------------------------------------------------------------------------------
// Test equivalence with a FieldBase.
//------------------------------------------------------------------------------
template<typename Value>
struct CrappyFieldCompareMethod {
  static bool compare(const ManagedVector<Value>& lhs, 
                      const ManagedVector<Value>& rhs) {
    if (&lhs == &rhs) return true;
    if (lhs.size() != rhs.size()) return false;
    for (size_t i = 0; i < lhs.size(); i++) if (lhs[i] != rhs[i]) return false; 
    return true;
    //return lhs == rhs;
  }
};

template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::operator==(const FieldBase<Dimension>& rhs) const {
  if (this->name() != rhs.name()) return false;
  if (this->nodeListPtr() != rhs.nodeListPtr()) return false;
  try {
    const Field<Dimension, DataType>* rhsPtr = dynamic_cast<const Field<Dimension, DataType>*>(&rhs);
    if (rhsPtr == 0) return false;
    return CrappyFieldCompareMethod<DataType>::compare(FieldViewType::mDataArray, rhsPtr->mDataArray);
  } catch (const std::bad_cast &) {
    return false;
  }
}

//------------------------------------------------------------------------------
// Element access by Node ID iterator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
DataType&
Field<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) {
  CHECK(itr.nodeListPtr() == this->nodeListPtr());
  CHECK(itr.nodeID() >= 0 && itr.nodeID() < numElements());
  return FieldViewType::mDataArray[itr.nodeID()];
}

template<typename Dimension, typename DataType>
inline
const DataType&
Field<Dimension, DataType>::
operator()(const NodeIteratorBase<Dimension>& itr) const {
  CHECK(itr.nodeListPtr() == this->nodeListPtr());
  CHECK(itr.nodeID() >= 0 && itr.nodeID() < numElements());
  return FieldViewType::mDataArray[itr.nodeID()];
}

template<typename Dimension, typename DataType>
inline
unsigned 
Field<Dimension, DataType>::numInternalElements() const {
  return this->nodeList().numInternalNodes();
}

template<typename Dimension, typename DataType>
inline
unsigned 
Field<Dimension, DataType>::numGhostElements() const {
  return this->nodeList().numGhostNodes();
}

//------------------------------------------------------------------------------
// Zero out the field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::Zero() {
  REQUIRE(mValid);
  std::fill(FieldViewType::begin(), FieldViewType::end(), DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// Addition with another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator+(const Field<Dimension, DataType>& rhs) const {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) += rhs(i);
  return result;
}

//------------------------------------------------------------------------------
// Subtract another field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator-(const Field<Dimension, DataType>& rhs) const {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) -= rhs(i);
  return result;
}

//------------------------------------------------------------------------------
// Add a single value to every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator+(const DataType& rhs) const {
  REQUIRE(valid());
  Field<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) += rhs;
  return result;
}

//------------------------------------------------------------------------------
// Subtract a single value from every element of a field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::operator-(const DataType& rhs) const {
  REQUIRE(valid());
  Field<Dimension, DataType> result(*this);
  const unsigned n = this->numElements();
  for (unsigned i = 0; i != n; ++i) result(i) -= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar Field
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator*(const Field<Dimension, Scalar>& rhs) const {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division by a Scalar Field
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator/(const Field<Dimension, Scalar>& rhs) const {
  REQUIRE(valid() && rhs.valid());
  REQUIRE(this->nodeListPtr() == rhs.nodeListPtr());
  Field<Dimension, DataType> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a Scalar
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator*(const Scalar& rhs) const {
  REQUIRE(valid());
  Field<Dimension, DataType> result(*this);
  result *= rhs;
  return result;
}

//------------------------------------------------------------------------------
// Division by a Scalar
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
Field<Dimension, DataType>
Field<Dimension, DataType>::
operator/(const Scalar& rhs) const {
  REQUIRE(valid());
  REQUIRE(rhs != 0.0);
  Field<Dimension, DataType> result(*this);
  result /= rhs;
  return result;
}

//------------------------------------------------------------------------------
// operator==(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator==(const Field<Dimension, DataType>& rhs) const {
  const auto n = this->size();
  if (n != rhs.size()) return false;
  auto result = true;
  auto i = 0;
  while (i < (int)n and result) {
    result = (*this)(i) == rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// operator!=(Field)
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
operator!=(const Field<Dimension, DataType>& rhs) const {
  return !((*this) == rhs);
}

//------------------------------------------------------------------------------
// Test if the field is valid.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::valid() const {
  //return true;
  return mValid && this->nodeListPtr();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::internalBegin() {
  return FieldViewType::mDataArray.begin();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::internalEnd() {
  CHECK(this->nodeList().firstGhostNode() >= 0 &&
         this->nodeList().firstGhostNode() <= this->nodeList().numNodes());
  return FieldViewType::mDataArray.begin() + this->nodeList().firstGhostNode();
}

//------------------------------------------------------------------------------
// Iterator pointing to the beginning of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::ghostBegin() {
  return this->internalEnd();
}

//------------------------------------------------------------------------------
// Iterator pointing to the end of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::iterator
Field<Dimension, DataType>::ghostEnd() {
  return FieldViewType::mDataArray.end();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::internalBegin() const {
  return FieldViewType::mDataArray.begin();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the internal field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::internalEnd() const {
  REQUIRE(this->nodeList().firstGhostNode() >= 0 &&
          this->nodeList().firstGhostNode() <= this->nodeList().numNodes());
  return FieldViewType::mDataArray.begin() + this->nodeList().firstGhostNode();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the beginning of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::ghostBegin() const {
  return this->internalEnd();
}

//------------------------------------------------------------------------------
// Const_iterator pointing to the end of the ghost field values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
typename Field<Dimension, DataType>::const_iterator
Field<Dimension, DataType>::ghostEnd() const {
  return FieldViewType::mDataArray.end();
}

//------------------------------------------------------------------------------
// Set the NodeList with which this field is defined.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::setNodeList(const NodeList<Dimension>& nodeList) {
  unsigned oldSize = this->size();
  this->setFieldBaseNodeList(nodeList);
  FieldViewType::mDataArray.resize(nodeList.numNodes());
  if (this->size() > oldSize) {
    for (unsigned i = oldSize; i < this->size(); ++i) {
      (*this)(i) = DataTypeTraits<DataType>::zero();
    }
  }
  mValid = true;
}

//------------------------------------------------------------------------------
// Resize the field to the given number of nodes.  This operation ignores
// the distinction between internal and ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeField(unsigned size) {
  REQUIRE(size == this->nodeList().numNodes());
  unsigned oldSize = this->size();
  FieldViewType::mDataArray.resize(size);
  if (oldSize < size) {
    std::fill(FieldViewType::mDataArray.begin() + oldSize,
              FieldViewType::mDataArray.end(),
              DataTypeTraits<DataType>::zero());
  }
  mValid = true;
}

//------------------------------------------------------------------------------
// Delete the given element id.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::deleteElement(int nodeID) {
  const unsigned originalSize = this->size();
  CONTRACT_VAR(originalSize);
  REQUIRE(nodeID >= 0 && nodeID < (int)originalSize);
  FieldViewType::mDataArray.erase(FieldViewType::mDataArray.begin() + nodeID);
  ENSURE(FieldViewType::mDataArray.size() == originalSize - 1);
}

//------------------------------------------------------------------------------
// Delete the given set of elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::deleteElements(const std::vector<int>& nodeIDs) {
  // The standalone method does the actual work.
  removeElements(FieldViewType::mDataArray, nodeIDs);
}

//------------------------------------------------------------------------------
// Serialize the chosen Field values onto a buffer
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<char>
Field<Dimension, DataType>::
packValues(const std::vector<int>& nodeIDs) const {
  return packFieldValues(*this, nodeIDs);
}

//------------------------------------------------------------------------------
// Unpack the given buffer into the requested field elements.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
unpackValues(const std::vector<int>& nodeIDs,
             const std::vector<char>& buffer) {
  unpackFieldValues(*this, nodeIDs, buffer);
}

//------------------------------------------------------------------------------
// Resize the field to the given number of internal nodes, preserving any ghost
// values at the end.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeFieldInternal(const unsigned size,
                                                const unsigned oldFirstGhostNode) {
  const unsigned currentSize = this->size();
  const unsigned currentInternalSize = oldFirstGhostNode;
  const unsigned numGhostNodes = this->nodeList().numGhostNodes();
  const unsigned newSize = size + numGhostNodes;
  REQUIRE(numGhostNodes == currentSize - oldFirstGhostNode);
  REQUIRE(newSize == this->nodeList().numNodes());

  // If there is ghost data, we must preserve it.
  std::vector<DataType,DataAllocator<DataType>> oldGhostValues(numGhostNodes);
  if (numGhostNodes > 0) {
    for (auto i = 0u; i != numGhostNodes; ++i) {
      const int j = oldFirstGhostNode + i;
      CHECK(i < numGhostNodes);
      CHECK(j < (int)this->size());
      oldGhostValues[i] = (*this)(j);
    }
  }

  // Resize the field data.
  FieldViewType::mDataArray.resize(newSize);

  // Fill in any new internal values.
  if (newSize > currentSize) {
    CHECK(currentInternalSize < this->nodeList().firstGhostNode());
    std::fill(FieldViewType::mDataArray.begin() + currentInternalSize,
              FieldViewType::mDataArray.begin() + this->nodeList().firstGhostNode(),
              DataTypeTraits<DataType>::zero());
  }

  // Fill the ghost data back in.
  if (numGhostNodes > 0) {
    for (auto i = 0u; i != numGhostNodes; ++i) {
      const int j = this->nodeList().firstGhostNode() + i;
      CHECK(i < oldGhostValues.size());
      CHECK(j < (int)this->size());
      (*this)(j) = oldGhostValues[i];
    }
  }

  mValid = true;
}

//------------------------------------------------------------------------------
// Resize the field to the given number of ghost nodes.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::resizeFieldGhost(const unsigned size) {
  const unsigned currentSize = this->size();
  const unsigned numInternalNodes = this->nodeList().numInternalNodes();
  const unsigned currentNumGhostNodes = currentSize - numInternalNodes;
  REQUIRE(currentNumGhostNodes >= 0);
  const unsigned newSize = numInternalNodes + size;
  REQUIRE(newSize == this->nodeList().numNodes());

  // Resize the field data.
  FieldViewType::mDataArray.resize(newSize);
  CHECK(this->size() == (newSize));

  // Fill in any new ghost values.
  if (newSize > currentSize) {
    std::fill(FieldViewType::mDataArray.begin() + numInternalNodes + currentNumGhostNodes,
              FieldViewType::mDataArray.end(),
              DataTypeTraits<DataType>::zero());
  }

  mValid = true;
}

//------------------------------------------------------------------------------
// Copy values between sets of indices.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
copyElements(const std::vector<int>& fromIndices,
             const std::vector<int>& toIndices) {
  REQUIRE(fromIndices.size() == toIndices.size());
  REQUIRE(std::all_of(fromIndices.begin(), fromIndices.end(),
                      [&](const int i) { return i >= 0 and i < (int)this->size(); }));
  REQUIRE(std::all_of(toIndices.begin(), toIndices.end(),
                      [&](const int i) { return i >= 0 and i < (int)this->size(); }));
  const auto ni = fromIndices.size();
  for (auto k = 0u; k < ni; ++k) (*this)(toIndices[k]) = (*this)(fromIndices[k]);
}

//------------------------------------------------------------------------------
// fixedSizeDataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
bool
Field<Dimension, DataType>::
fixedSizeDataType() const {
  return DataTypeTraits<DataType>::fixedSize();
}

//------------------------------------------------------------------------------
// numValsInDataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
int
Field<Dimension, DataType>::
numValsInDataType() const {
  return DataTypeTraits<DataType>::numElements(DataType());
}

//------------------------------------------------------------------------------
// sizeofDataType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
int
Field<Dimension, DataType>::
sizeofDataType() const {
  return sizeof(DataTypeTraits<DataType>::zero());
}

//------------------------------------------------------------------------------
// computeCommBufferSize
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
int
Field<Dimension, DataType>::
computeCommBufferSize(const std::vector<int>& packIndices,
                      const int sendProc,
                      const int recvProc) const {
  return computeBufferSize(*this, packIndices, sendProc, recvProc);
}

//------------------------------------------------------------------------------
// Serialize the Field into a vector<char>.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<char>
Field<Dimension, DataType>::
serialize() const {
  const size_t n = numInternalElements();
  vector<char> buf;
  packElement(this->name(), buf);
  packElement(n, buf);
  for (auto i = 0u; i < n; ++i) packElement((*this)[i], buf);
  return buf;
}

//------------------------------------------------------------------------------
// Deserialize the values from a vector<char> into this Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
void
Field<Dimension, DataType>::
deserialize(const std::vector<char>& buf) {
  auto itr = buf.begin();
  std::string nm;
  unpackElement(nm, itr, buf.end());
  this->name(nm);
  size_t n;
  unpackElement(n, itr, buf.end());
  VERIFY2(n == this->numInternalElements(),
          "Field ERROR: attempt to deserialize wrong number of elements: " << n << " != " << this->numInternalElements());
  for (auto i = 0u; i < n; ++i) unpackElement((*this)[i], itr, buf.end());
}

//------------------------------------------------------------------------------
// Construct std::vectors of pointers to the values.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
internalValues() const {
  std::vector<DataType> result;
  result.reserve(this->nodeList().numInternalNodes());
  for (const_iterator itr = internalBegin();
       itr != internalEnd();
       ++itr) result.push_back(*itr);
  ENSURE(result.size() == this->nodeList().numInternalNodes());
  return result;
}

template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
ghostValues() const {
  std::vector<DataType> result;
  result.reserve(this->nodeList().numGhostNodes());
  for (const_iterator itr = ghostBegin();
       itr != ghostEnd();
       ++itr) result.push_back(*itr);
  ENSURE(result.size() == this->nodeList().numGhostNodes());
  return result;
}

template<typename Dimension, typename DataType>
inline
std::vector<DataType>
Field<Dimension, DataType>::
allValues() const {
  std::vector<DataType> result;
  result.reserve(this->nodeList().numNodes());
  for (const_iterator itr = FieldViewType::begin();
       itr != FieldViewType::end();
       ++itr) result.push_back(*itr);
  ENSURE(result.size() == this->nodeList().numNodes());
  return result;
}

//****************************** Global Functions ******************************
//------------------------------------------------------------------------------
// Multiplication by another Field.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const Field<Dimension, DataType>& lhs,
          const Field<Dimension, OtherDataType>& rhs) {
  CHECK(lhs.valid() && rhs.valid());
  CHECK(lhs.nodeList().numNodes() == rhs.nodeList().numNodes());
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, DataType>&>(lhs).nodeList());
  for (auto i = 0u; i < result.numElements(); ++i) {
    result(i) = lhs(i) * rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Multiplication by a single value.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const Field<Dimension, DataType>& lhs,
          const OtherDataType& rhs) {
  CHECK(lhs.valid());
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, DataType>&>(lhs).nodeList());
  for (int i = 0; i < result.numElements(); ++i) {
    result(i) = lhs(i) * rhs;
  }
  return result;
}

template<typename Dimension, typename DataType, typename OtherDataType>
Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
operator*(const DataType& lhs,
          const Field<Dimension, OtherDataType>& rhs) {
  CHECK(rhs.valid());
  Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
    result("product", const_cast<Field<Dimension, OtherDataType>&>(rhs).nodeList());
  for (auto i = 0u; i < result.numElements(); ++i) {
    result(i) = lhs * rhs(i);
  }
  return result;
}

//------------------------------------------------------------------------------
// Input (istream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::istream&
operator>>(std::istream& is, Field<Dimension, DataType>& field) {

  // Start by reading the number of elements.
  int numElementsInStream;
  is >> numElementsInStream;
  CHECK(numElementsInStream == (int)field.nodeList().numInternalNodes());

  // Read in the elements.
  for (typename Field<Dimension, DataType>::iterator itr = field.internalBegin();
       itr < field.internalEnd();
       ++itr) {
    is >> *itr;
  }
  return is;
}

//------------------------------------------------------------------------------
// Output (ostream) operator.
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
std::ostream&
operator<<(std::ostream& os, const Field<Dimension, DataType>& field) {

  // Write the number of internal elements.
  os << field.nodeList().numInternalNodes() << " ";

  // Write the internal elements.
  for (typename Field<Dimension, DataType>::const_iterator itr = field.internalBegin();
       itr < field.internalEnd();
       ++itr) {
    os << *itr << " ";
  }
//   os << endl;
  return os;
}

//------------------------------------------------------------------------------
// getAxomType
//------------------------------------------------------------------------------
template<typename Dimension, typename DataType>
inline
axom::sidre::DataTypeId
Field<Dimension, DataType>::
getAxomTypeID() const {
  return DataTypeTraits<DataType>::axomTypeID();
}

} // namespace Spheral
