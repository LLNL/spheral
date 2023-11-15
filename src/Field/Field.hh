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

  //using ContainerType = ManagedVector<DataType>;
  using FieldViewType = FieldView<Dimension, DataType>;
  using ContainerType = typename FieldViewType::ContainerType;

  using iterator = typename ContainerType::iterator;
  using const_iterator = typename ContainerType::const_iterator;
  //typedef typename std::vector<DataType,DataAllocator<DataType>>::iterator iterator;
  //typedef typename std::vector<DataType,DataAllocator<DataType>>::const_iterator const_iterator;

  // Constructors.
  SPHERAL_HOST explicit Field(FieldName name);
  SPHERAL_HOST Field(FieldName name, const Field& field);
  SPHERAL_HOST Field(FieldName name,
                     const NodeList<Dimension>& nodeList);
  SPHERAL_HOST Field(FieldName name,
                     const NodeList<Dimension>& nodeList,
                     DataType value);
  SPHERAL_HOST Field(FieldName name,
                     const NodeList<Dimension>& nodeList, 
                     const ContainerType& array);
  SPHERAL_HOST Field(const NodeList<Dimension>& nodeList, const Field& field);
  
  SPHERAL_HOST Field(const Field& field);
  SPHERAL_HOST virtual std::shared_ptr<FieldBase<Dimension> > clone() const override;

  // Destructor.
  SPHERAL_HOST ~Field();

  // Assignment operator.
  SPHERAL_HOST virtual FieldBase<Dimension>& operator=(const FieldBase<Dimension>& rhs) override;
  SPHERAL_HOST Field& operator=(const Field& rhs);
  SPHERAL_HOST Field& operator=(const ContainerType& rhs);
  SPHERAL_HOST Field& operator=(const DataType& rhs);

  SPHERAL_HOST inline FieldViewType toView() const { return FieldViewType(*this); } 

  // Required method to test equivalence with a FieldBase.
  SPHERAL_HOST virtual bool operator==(const FieldBase<Dimension>& rhs) const override;

  SPHERAL_HOST DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  SPHERAL_HOST const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  // The number of elements in the field.
  SPHERAL_HOST unsigned numInternalElements() const;
  SPHERAL_HOST unsigned numGhostElements() const;

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
  SPHERAL_HOST bool valid() const;

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

  // **************************************************************************
  // *** Explicit Forward Interface From FieldView ***
  // **************************************************************************

  using FieldViewType::operator();
  using FieldViewType::operator==;
  using FieldViewType::operator!=;
  
  //// Element access.
  //DataType& operator[](const unsigned int index) { return FieldViewType::operator[](index); }
  //DataType& operator[](const unsigned int index) const { return FieldViewType::operator[](index); }

  //SPHERAL_HOST inline DataType& operator()(int index) { return FieldViewType::operator()(index); }
  //SPHERAL_HOST inline const DataType& operator()(int index) const { return FieldViewType::operator()(index); }
  ////using FieldViewType::operator();
  //SPHERAL_HOST inline DataType& at(int index) { return FieldViewType::at(index); }
  //SPHERAL_HOST inline const DataType& at(int index) const { return FieldViewType::at(index); }
  ////using FieldViewType::at;
  //SPHERAL_HOST inline unsigned numElements() const { return FieldViewType::numElements(); }
  ////using FieldViewType::numElements;

  //// Methods to apply limits to Field data members.
  //inline void applyMin(const DataType& dataMin) { FieldViewType::applyMin(dataMin); }
  //inline void applyMax(const DataType& dataMax) { FieldViewType::applyMax(dataMax); }

  //inline void applyScalarMin(const double dataMin) { FieldViewType::applyScalarMin(dataMin); }
  //inline void applyScalarMax(const double dataMax) { FieldViewType::applyScalarMax(dataMax); }

  //// Some useful reduction operations.
  //inline DataType sumElements() const { return FieldViewType::sumElements(); }
  //inline DataType min() const { return FieldViewType::min(); }
  //inline DataType max() const { return FieldViewType::max(); }

  //// Some useful reduction operations (local versions -- no MPI reductions)
  //inline DataType localSumElements() const { return FieldViewType::localSumElements(); }
  //inline DataType localMin() const { return FieldViewType::localMin(); }
  //inline DataType localMax() const { return FieldViewType::localMax(); }

  //bool operator>(const Field& rhs) const { return FieldViewType::operator>(rhs); }
  //bool operator<(const Field& rhs) const { return FieldViewType::operator<(rhs); }
  //bool operator>=(const Field& rhs) const { return FieldViewType::operator>=(rhs); }
  //bool operator<=(const Field& rhs) const { return FieldViewType::operator<=(rhs); }

  //// Comparison operators (Field-value element wise).
  //bool operator==(const DataType& rhs) const { return FieldViewType::operator==(rhs); }
  //bool operator!=(const DataType& rhs) const { return FieldViewType::operator!=(rhs); }
  //bool operator>(const DataType& rhs) const { return FieldViewType::operator>(rhs); }
  //bool operator<(const DataType& rhs) const { return FieldViewType::operator<(rhs); }
  //bool operator>=(const DataType& rhs) const { return FieldViewType::operator>=(rhs); }
  //bool operator<=(const DataType& rhs) const { return FieldViewType::operator<=(rhs); }

  //// Provide the standard iterator methods over the field.
  //iterator begin() { return FieldViewType::begin(); }
  //iterator end() { return FieldViewType::end(); }

  //const_iterator begin() const { return FieldViewType::begin(); }
  //const_iterator end() const { return FieldViewType::end(); }

  //SPHERAL_HOST virtual unsigned size() const override { return FieldViewType::size(); }

private:
  //--------------------------- Private Interface ---------------------------//
  // No default constructor.
  Field();
  // Private Data
  bool mValid;

  //using FieldViewType::operatorField();
  //operator Field<Dimension, DataType>() {return Field<Dimension, DataType>(*this); }
};

//template<typename Dimension, typename DataType>
//Field<Dimension, DataType> deepCopy(Field<Dimension, DataType> const& rhs)

} // namespace Spheral

#include "FieldInline.hh"

#else

// Forward declare the Field class.
namespace Spheral {
  template<typename Dimension, typename DataType> class Field;
} // namespace Spheral

#endif
