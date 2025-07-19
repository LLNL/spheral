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

#include "Utilities/Logger.hh"

#include "axom/sidre.hpp"
#include "chai/ExecutionSpaces.hpp"
#include "chai/ManagedArray.hpp"
#include "chai/PointerRecord.hpp"
#include "chai/Types.hpp"

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
    public FieldBase<Dimension> {
   
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

  typedef typename std::vector<DataType,DataAllocator<DataType>>::iterator iterator;
  typedef typename std::vector<DataType,DataAllocator<DataType>>::const_iterator const_iterator;

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
  virtual ~Field();

  // Assignment operator.
  virtual FieldBase<Dimension>& operator=(const FieldBase<Dimension>& rhs) override;
  Field& operator=(const Field& rhs);
  Field& operator=(const std::vector<DataType,DataAllocator<DataType>>& rhs);
  Field& operator=(const DataType& rhs);

  // Required method to test equivalence with a FieldBase.
  virtual bool operator==(const FieldBase<Dimension>& rhs) const override;

  // Element access.
  DataType& operator()(int index);
  const DataType& operator()(int index) const;

  DataType& operator()(const NodeIteratorBase<Dimension>& itr);
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const;

  DataType& at(int index);
  const DataType& at(int index) const;

  // The number of elements in the field.
  unsigned numElements() const;
  unsigned numInternalElements() const;
  unsigned numGhostElements() const;
  virtual unsigned size() const override;

  // Zero out the field elements.
  virtual void Zero() override;

  // Methods to apply limits to Field data members.
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const double dataMin);
  void applyScalarMax(const double dataMax);

  // Standard field additive operators.
  Field operator+(const Field& rhs) const;
  Field operator-(const Field& rhs) const;

  Field& operator+=(const Field& rhs);
  Field& operator-=(const Field& rhs);

  Field operator+(const DataType& rhs) const;
  Field operator-(const DataType& rhs) const;

  Field& operator+=(const DataType& rhs);
  Field& operator-=(const DataType& rhs);

//   // Multiplication of two fields, possibly by another DataType.
//   template<typename OtherDataType>
//   Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
//   operator*(const Field<Dimension, OtherDataType>& rhs) const;

//   template<typename OtherDataType>
//   Field<Dimension, typename CombineTypes<DataType, OtherDataType>::ProductType>
//   operator*(const OtherDataType& rhs) const;

  // Multiplication and division by scalar(s)
  Field<Dimension, DataType> operator*(const Field<Dimension, Scalar>& rhs) const;
  Field<Dimension, DataType> operator/(const Field<Dimension, Scalar>& rhs) const;

  Field<Dimension, DataType>& operator*=(const Field<Dimension, Scalar>& rhs);
  Field<Dimension, DataType>& operator/=(const Field<Dimension, Scalar>& rhs);

  Field<Dimension, DataType> operator*(const Scalar& rhs) const;
  Field<Dimension, DataType> operator/(const Scalar& rhs) const;

  Field<Dimension, DataType>& operator*=(const Scalar& rhs);
  Field<Dimension, DataType>& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const;
  DataType localMin() const;
  DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  bool operator==(const Field& rhs) const;
  bool operator!=(const Field& rhs) const;
  bool operator>(const Field& rhs) const;
  bool operator<(const Field& rhs) const;
  bool operator>=(const Field& rhs) const;
  bool operator<=(const Field& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator>(const DataType& rhs) const;
  bool operator<(const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

//   // Interpolate from this Field onto the given position.  Assumes that the
//   // neighbor initializations have already been performed for the given
//   // position!
//   DataType operator()(const Vector& r,
//                       const TableKernel<Dimension>& W) const;

//   // Interpolate from this Field onto a new Field defined at the positions
//   // of the given NodeList.
//   Field<Dimension, DataType>
//   sampleField(const NodeList<Dimension>& splatNodeList,
//               const TableKernel<Dimension>& W) const;

//   // Conservatively splat values from this Field onto a new Field defined
//   // at the positions of the given NodeList, using the MASH formalism.
//   Field<Dimension, DataType>
//   splatToFieldMash(const NodeList<Dimension>& splatNodeList,
//                    const TableKernel<Dimension>& W) const;

  // Test if this Field is in a valid, internally consistent state.
  bool valid() const;

  // Provide the standard iterator methods over the field.
  iterator begin();
  iterator end();
  iterator internalBegin();
  iterator internalEnd();
  iterator ghostBegin();
  iterator ghostEnd();

  const_iterator begin() const;
  const_iterator end() const;
  const_iterator internalBegin() const;
  const_iterator internalEnd() const;
  const_iterator ghostBegin() const;
  const_iterator ghostEnd() const;

  // Index operator.
  DataType& operator[](const unsigned int index);
  const DataType& operator[](const unsigned int index) const;

  // Required functions from FieldBase
  virtual void setNodeList(const NodeList<Dimension>& nodeList) override;
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

  using ViewType = FieldView<Dimension, DataType>;

  ViewType toView()
  {
    auto func = [](
        const chai::PointerRecord *,
        chai::Action,
        chai::ExecutionSpace) {};

    return this->toView(func);
  }

  template<typename T=DataType, typename F>
  std::enable_if_t<std::is_trivially_copyable<T>::value, ViewType>
  toView(F&& extension)
  {
    if (mManagedData.size() != mDataArray.size() ||
        mManagedData.data(chai::CPU, false) != mDataArray.data()) {

      mManagedData.free();

      mManagedData = chai::makeManagedArray(
          mDataArray.data(), mDataArray.size(), chai::CPU, false);

      mManagedData.setUserCallback(
        getFieldCallback(std::forward<F>(extension)));
    }
    return ViewType(mManagedData);
  }

  template<typename T=DataType, typename F>
  std::enable_if_t<!std::is_trivially_copyable<T>::value, ViewType>
  toView(F&&)
  {
    ASSERT2(false, "Spheral::Field::toView() Is invalid when Field::DataType is not trivially copyable.");
    return ViewType(mManagedData);
  }

protected:
  virtual void resizeField(unsigned size) override;
  virtual void resizeFieldInternal(unsigned size, unsigned oldFirstGhostNode) override;
  virtual void resizeFieldGhost(unsigned size) override;
  virtual void deleteElement(int nodeID) override;
  virtual void deleteElements(const std::vector<int>& nodeIDs) override;

  template<typename F>
  auto getFieldCallback(F callback) {
    return [n = this->name(), callback](
      const chai::PointerRecord * record,
      chai::Action action,
      chai::ExecutionSpace space) {
        if (action == chai::ACTION_MOVE) {
          if (space == chai::CPU)
            DEBUG_LOG << "Field :" << n << ": MOVED to the CPU";
          if (space == chai::GPU)
            DEBUG_LOG << "Field :" << n << ": MOVED to the GPU";
        }
        else if (action == chai::ACTION_ALLOC) {
          if (space == chai::CPU)
            DEBUG_LOG << "Field :" << n << ": ALLOC on the CPU";
          if (space == chai::GPU)
            DEBUG_LOG << "Field :" << n << ": ALLOC on the GPU";
        }
        else if (action == chai::ACTION_FREE) {
          if (space == chai::CPU)
            DEBUG_LOG << "Field :" << n << ": FREE on the CPU";
          if (space == chai::GPU)
            DEBUG_LOG << "Field :" << n << ": FREE on the GPU";
        }
        callback(record, action, space);
      };
  };

private:
  //--------------------------- Private Interface ---------------------------//
  // Private Data
  std::vector<DataType, DataAllocator<DataType>> mDataArray;

  // ManagedArray owned by the field to ensure proper lifetime of the GPU data.
  chai::ManagedArray<DataType> mManagedData;


  bool mValid;

  // No default constructor.
  Field();
};

} // namespace Spheral

#include "FieldInline.hh"

#endif
