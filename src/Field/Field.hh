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
#ifndef __Spheral_Field_hh__
#define __Spheral_Field_hh__

#include "FieldBase.hh"
#include "FieldView.hh"
#include "axom/sidre.hpp"
#include "SphArray.hh"
#include "config.hh"

#include <vector>

#ifdef USE_UVM
#include "uvm_allocator.hh"
#endif

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
    public FieldBase<Dimension>, public FieldView<Dimension, DataType> {

public:
  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  
  typedef typename FieldBase<Dimension>::FieldName FieldName;
  typedef Dimension FieldDimension;
  typedef DataType FieldDataType;
  typedef DataType value_type;      // STL compatibility.

  using FieldViewType = FieldView<Dimension, DataType>;
  using ContainerType = typename FieldViewType::ContainerType;

  using iterator = typename ContainerType::iterator;
  using const_iterator = typename ContainerType::const_iterator;

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
                     const ContainerType& array);
  Field(const NodeList<Dimension>& nodeList, const Field& field);
  
  Field(const Field& field);
  virtual std::shared_ptr<FieldBase<Dimension> > clone() const override;

  // Destructor.
  ~Field();

  // Assignment operator.
  virtual FieldBase<Dimension>& operator=(const FieldBase<Dimension>& rhs) override;
  Field& operator=(const Field& rhs);
  Field& operator=(const ContainerType& rhs);
  Field& operator=(const DataType& rhs);

  inline FieldViewType toView() const { return FieldViewType(*this); } 

  // Required method to test equivalence with a FieldBase.
  virtual bool operator==(const FieldBase<Dimension>& rhs) const override;

  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  // The number of elements in the field.
  unsigned numInternalElements() const;
  unsigned numGhostElements() const;

  // Standard field additive operators.
  Field operator+(const Field& rhs) const;
  Field operator-(const Field& rhs) const;

  inline Field& operator+=(const Field& rhs) { FieldViewType::operator+=(rhs); return *this; }
  inline Field& operator-=(const Field& rhs) { FieldViewType::operator-=(rhs); return *this; }

  Field operator+(const DataType& rhs) const;
  Field operator-(const DataType& rhs) const;

  inline Field& operator+=(const DataType& rhs) { FieldViewType::operator+=(rhs); return *this; }
  inline Field& operator-=(const DataType& rhs) { FieldViewType::operator-=(rhs); return *this; }

  // Multiplication and division by scalar(s)
  Field operator*(const Field<Dimension, Scalar>& rhs) const;
  Field operator/(const Field<Dimension, Scalar>& rhs) const;

  inline Field& operator*=(const Field<Dimension, Scalar>& rhs) { FieldViewType::operator*=(rhs); return *this; }
  inline Field& operator/=(const Field<Dimension, Scalar>& rhs) { FieldViewType::operator/=(rhs); return *this; }

  Field operator*(const Scalar& rhs) const;
  Field operator/(const Scalar& rhs) const;

  inline Field& operator*=(const Scalar& rhs) { FieldViewType::operator*=(rhs); return *this; }
  inline Field& operator/=(const Scalar& rhs) { FieldViewType::operator/=(rhs); return *this; }

  // Comparison operators (Field-Field element wise).
  bool operator==(const Field& rhs) const;
  bool operator!=(const Field& rhs) const;

  iterator internalBegin();
  iterator internalEnd();
  iterator ghostBegin();
  iterator ghostEnd();

  const_iterator internalBegin() const;
  const_iterator internalEnd() const;
  const_iterator ghostBegin() const;
  const_iterator ghostEnd() const;

  // Test if this Field is in a valid, internally consistent state.
  bool valid() const;

  // Required functions from FieldBase
  virtual void Zero() override;
  virtual void setNodeList(const NodeList<Dimension>& nodeList) override;
  virtual void resizeField(unsigned size) override;
  virtual void resizeFieldInternal(unsigned size, unsigned oldFirstGhostNode) override;
  virtual void resizeFieldGhost(unsigned size) override;
  virtual void deleteElement(int nodeID) override;
  virtual void deleteElements(const std::vector<int>& nodeIDs) override;
  virtual std::vector<char> packValues(const std::vector<int>& nodeIDs) const override;
  virtual void unpackValues(const std::vector<int>& nodeIDs,
                            const std::vector<char>& buffer) override;
  virtual void copyElements(const std::vector<int>& fromIndices,
                            const std::vector<int>& toIndices) override;
  virtual bool fixedSizeDataType() const override;
  virtual int numValsInDataType() const override;
  virtual int sizeofDataType() const override;
  virtual int computeCommBufferSize(const std::vector<int>& packIndices,
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

  // Explicit Forward Interface From FieldView for overloaded methods. 
  using FieldViewType::operator();
  using FieldViewType::operator==;
  using FieldViewType::operator!=;
  using FieldViewType::numElements;
  
private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor.
  Field();
  // Private Data
  bool mValid;
};

} // namespace Spheral

#include "FieldInline.hh"

#else

// Forward declare the Field class.
namespace Spheral {
  template<typename Dimension, typename DataType> class Field;
} // namespace Spheral

#endif
