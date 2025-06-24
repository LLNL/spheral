//---------------------------------Spheral++----------------------------------//
// Field -- Provide a field of a type (Scalar, Vector, Tensor) over the nodes
//          in a NodeList.
//
// This version of the Field is based on standard constructs like the STL
// vector.  This will certainly be slower at run time than the Blitz Array
// class.
//
// Created by JMO, Thu Jun 10 23:26:50 PDT 1999
//----------------------------------------------------------------------------//
#ifndef __Spheral_Field__
#define __Spheral_Field__

#include "Field/FieldBase.hh"
#include "Field/FieldSpan.hh"

#include "axom/sidre.hpp"

#ifdef USE_UVM
#include "uvm_allocator.hh"
#endif

#include <vector>

namespace Spheral {

template<typename Dimension> class NodeIteratorBase;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class NodeList;
template<typename Dimension> class TableKernel;

#ifdef USE_UVM
template<typename DataType>
using DataAllocator = typename uvm_allocator::UVMAllocator<DataType>;
#else
template<typename DataType>
using DataAllocator = std::allocator<DataType>;
#endif

template<typename Dimension, typename DataType>
class Field:
    public FieldBase<Dimension>,
    public FieldSpan<Dimension, DataType> {
   
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  using Vector = typename Dimension::Vector;
  using Tensor = typename Dimension::Tensor;
  using SymTensor = typename Dimension::SymTensor;
  
  using FieldName = typename FieldBase<Dimension>::FieldName;
  using FieldDimension = Dimension;
  using FieldDataType = DataType;
  using value_type = DataType;      // STL compatibility.

  using iterator = typename FieldSpan<Dimension, DataType>::iterator;
  using const_iterator = typename std::vector<DataType,DataAllocator<DataType>>::const_iterator;

  using ViewType = FieldSpan<Dimension, DataType>;

  // Bring in various methods hidden in FieldSpan
  using FieldSpan<Dimension, DataType>::operator();
  using FieldSpan<Dimension, DataType>::operator[];

  // Constructors.
  explicit Field(FieldName name);
  Field(FieldName name, const Field& field);
  Field(FieldName name,
        const NodeList<Dimension>& nodeList);
  Field(FieldName name,
        const NodeList<Dimension>& nodeList,
        DataType value);
  Field(FieldName name,
        const NodeList<Dimension>& nodeList, 
        const std::vector<DataType,DataAllocator<DataType>>& array);
  Field(const NodeList<Dimension>& nodeList, const Field& field);
  Field(const Field& field);
  virtual std::shared_ptr<FieldBase<Dimension> > clone() const override;

  // Destructor.
  virtual ~Field() = default;

  // Assignment operator.
  virtual FieldBase<Dimension>& operator=(const FieldBase<Dimension>& rhs) override;
  Field& operator=(const Field& rhs);
  Field& operator=(const std::vector<DataType,DataAllocator<DataType>>& rhs);
  Field& operator=(const DataType& rhs);

  // Comparisons
  virtual bool operator==(const FieldBase<Dimension>& rhs) const override;

  // Foward the FieldSpan Field-Field comparison operators
  bool operator==(const Field& rhs) const { return FieldSpan<Dimension, DataType>::operator==(rhs); }
  bool operator!=(const Field& rhs) const { return FieldSpan<Dimension, DataType>::operator!=(rhs); }

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const { return FieldSpan<Dimension, DataType>::operator==(rhs); }
  bool operator!=(const DataType& rhs) const { return FieldSpan<Dimension, DataType>::operator!=(rhs); }
  bool operator> (const DataType& rhs) const { return FieldSpan<Dimension, DataType>::operator> (rhs); }
  bool operator< (const DataType& rhs) const { return FieldSpan<Dimension, DataType>::operator< (rhs); }
  bool operator>=(const DataType& rhs) const { return FieldSpan<Dimension, DataType>::operator>=(rhs); }
  bool operator<=(const DataType& rhs) const { return FieldSpan<Dimension, DataType>::operator<=(rhs); }

  // Element access by NodeIterator
  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  // The number of elements in the field.
  virtual size_t size() const override                                      { return mDataArray.size(); }

  // Zero out the field elements.
  virtual void Zero() override;

  // Standard field additive operators.
  Field operator+(const Field& rhs) const;
  Field operator-(const Field& rhs) const;

  Field operator+(const DataType& rhs) const;
  Field operator-(const DataType& rhs) const;

  // Multiplication and division by scalar(s)
  Field operator*(const Field<Dimension, Scalar>& rhs) const;
  Field operator/(const Field<Dimension, Scalar>& rhs) const;

  Field operator*(const Scalar& rhs) const;
  Field operator/(const Scalar& rhs) const;

  // Provide the standard iterator methods over the field.
  const_iterator begin() const                                              { return mDataArray.begin(); }
  const_iterator end() const                                                { return mDataArray.end(); }
  const_iterator internalBegin() const                                      { return mDataArray.begin(); }
  const_iterator internalEnd() const                                        { return mDataArray.begin() + mNumInternalElements; }
  const_iterator ghostBegin() const                                         { return mDataArray.begin() + mNumInternalElements; }
  const_iterator ghostEnd() const                                           { return mDataArray.end(); }

  // We have to explicitly redefine the non-const iterators
  iterator begin()                                                          { return FieldSpan<Dimension, DataType>::begin(); }
  iterator end()                                                            { return FieldSpan<Dimension, DataType>::end(); }
  iterator internalBegin()                                                  { return FieldSpan<Dimension, DataType>::internalBegin(); }
  iterator internalEnd()                                                    { return FieldSpan<Dimension, DataType>::internalEnd(); }
  iterator ghostBegin()                                                     { return FieldSpan<Dimension, DataType>::ghostBegin(); }
  iterator ghostEnd()                                                       { return FieldSpan<Dimension, DataType>::ghostEnd(); }

  // Required functions from FieldBase
  virtual void setNodeList(const NodeList<Dimension>& nodeList) override;
  virtual std::vector<char> packValues(const std::vector<size_t>& nodeIDs) const override;
  virtual void unpackValues(const std::vector<size_t>& nodeIDs,
                            const std::vector<char>& buffer) override;
  virtual void copyElements(const std::vector<size_t>& fromIndices,
                            const std::vector<size_t>& toIndices) override;
  virtual bool fixedSizeDataType() const override;
  virtual size_t numValsInDataType() const override;
  virtual size_t sizeofDataType() const override;
  virtual size_t computeCommBufferSize(const std::vector<size_t>& packIndices,
                                       const int sendProc,
                                       const int recvProc) const override;

  // Serialization methods
  std::vector<char> serialize() const;
  void deserialize(const std::vector<char>& buf);

  // Provide std::vector copies of the data.  This is mostly useful for the
  // python interface.
  std::vector<DataType> internalValues() const;
  std::vector<DataType> ghostValues() const;
  std::vector<DataType> allValues() const;

  // Functions to help with storing the field in a Sidre datastore.
  axom::sidre::DataTypeId getAxomTypeID() const;

  // Return a view of the Field (appropriate for on accelerator devices)
  ViewType& view();

  // No default constructor.
  Field() = delete;

protected:
  virtual void resizeField(size_t size) override;
  virtual void resizeFieldInternal(size_t size, size_t oldFirstGhostNode) override;
  virtual void resizeFieldGhost(size_t size) override;
  virtual void deleteElement(size_t nodeID) override;
  virtual void deleteElements(const std::vector<size_t>& nodeIDs) override;

private:
  //--------------------------- Private Interface ---------------------------//
  // Private Data
  std::vector<DataType, DataAllocator<DataType>> mDataArray;

  friend FieldSpan<Dimension, DataType>;
  using FieldSpan<Dimension, DataType>::mDataSpan;
  using FieldSpan<Dimension, DataType>::mNumInternalElements;
  using FieldSpan<Dimension, DataType>::mNumGhostElements;
};

} // namespace Spheral

#include "FieldInline.hh"

#endif
