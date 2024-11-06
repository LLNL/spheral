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
#include "axom/sidre.hpp"
#include "Utilities/ValueViewInterface.hh"

#include <vector>

#ifdef USE_UVM
#include "uvm_allocator.hh"
#endif

namespace Spheral {

template<typename Dimension> class NodeIteratorBase;
template<typename Dimension> class CoarseNodeIterator;
template<typename Dimension> class RefineNodeIterator;
template<typename Dimension> class TableKernel;
template<typename Dimension> class NodeList;

#ifdef USE_UVM
template<typename DataType>
using DataAllocator = typename uvm_allocator::UVMAllocator<DataType>;
#else
template<typename DataType>
using DataAllocator = std::allocator<DataType>;
#endif

VVI_IMPL_BEGIN

template<typename Dimension, typename DataType>
class Field: 
    public FieldBase<Dimension> {
  //--------------------------- Private Interface ---------------------------//
  // Private Data
  using container_type = vvi::vector<DataType>; 
  container_type mDataArray;
  bool mValid;

  // No default constructor.
public:
  Field();
   
  VVI_IMPL_DEEPCOPY(Field, mDataArray)
  VVI_IMPL_COMPARE(Field, mDataArray)

  //--------------------------- Public Interface ---------------------------//
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  
  using FieldName = typename FieldBase<Dimension>::FieldName;
  typedef Dimension FieldDimension;
  typedef DataType FieldDataType;
  typedef DataType value_type;      // STL compatibility.

  typedef typename container_type::iterator iterator;
  typedef typename container_type::const_iterator const_iterator;

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
  virtual void setNodeList(const Spheral::NodeList<Dimension>& nodeList) override;
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

};

VVI_IMPL_END

#ifdef VVI_ENABLED

template<typename Dimension, typename DataType>
class Field;

#define FieldView__(code) PTR_VIEW_METACLASS_DECL((Field<Dimension, DataType>), (FieldView), (vvimpl::Field<Dimension, DataType>), (code))
#define Field__(code) PTR_VALUE_METACLASS_DECL((Field), (FieldView<Dimension, DataType>), (UNPACK code))


template<typename Dimension, typename DataType>
class FieldView__(
  // Field inherits from FieldBase so we would like to be able to implicitly
  // upcast a fieldview object to a fieldbaseview.
  VVI_UPCAST_CONVERSION_OP((FieldBaseView<Dimension>))

  DataType& operator()(int index) { return VVI_IMPL_INST().operator()(index); }
  const DataType& operator()(int index) const { return VVI_IMPL_INST().operator()(index); }

  DataType& operator()(const NodeIteratorBase<Dimension>& itr) { return VVI_IMPL_INST().operator()(itr); }
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const { return VVI_IMPL_INST().operator()(itr); }
);


template<typename Dimension, typename DataType>
class Field__((
public:
  using Scalar = typename ImplType::Scalar;
  using Vector = typename ImplType::Vector;
  using Tensor = typename ImplType::Tensor;
  using SymTensor = typename ImplType::SymTensor;
  
  using FieldName = typename ImplType::FieldName;
  using FieldDimension = typename ImplType::FieldDimension;
  using FieldDataType = typename ImplType::FieldDataType;
  using value_type = typename ImplType::value_type;

  using iterator = typename ImplType::iterator;
  using const_iterator = typename ImplType::const_iterator;

//  // Custom Ctor, note we need to create the underlying implementation 
//  // object on ctor of value interfaces.
  VVI_VALUE_DEF_CTOR(Field)
  Field(FieldName name) : //VVI_VALUE_CTOR_ARGS(name) {}
    Base(VVI_MAKE_SHARED<ImplType>(name)) {}

  Field(FieldName name,
        const NodeList<Dimension>& nodeList) :
    Base(VVI_MAKE_SHARED<ImplType>(name, nodeList)) {}

  Field(FieldName name,
        const NodeList<Dimension>& nodeList,
        DataType value) :
    Base(VVI_MAKE_SHARED<ImplType>(name, nodeList, value)) {}

//
//  // Value semantics dictate that we free underlying data upon destruction.
//  ~Field() { VVI_IMPL_INST()().mDataArray.free(); }

  // Assignment operator.
  //FieldBase<Dimension>& operator=(const FieldBase<Dimension>& rhs) override;
  Field& operator=(const Field& rhs) = default;
  //Field& operator=(const std::vector<DataType,DataAllocator<DataType>>& rhs);
  Field& operator=(const DataType& rhs) { return VVI_IMPL_INST().operator=(rhs); }

  // Required method to test equivalence with a FieldBase.
  //virtual bool operator==(const FieldBase<Dimension>& rhs) const override;

  // Element access.
  DataType& operator()(int index) { return VVI_IMPL_INST().operator()(index); }
  const DataType& operator()(int index) const { return VVI_IMPL_INST().operator()(index); }

  DataType& operator()(const NodeIteratorBase<Dimension>& itr) { return VVI_IMPL_INST().operator()(itr); }
  const DataType& operator()(const NodeIteratorBase<Dimension>& itr) const { return VVI_IMPL_INST().operator()(itr); }

  DataType& at(int index) { return VVI_IMPL_INST().at(index); }
  const DataType& at(int index) const { return VVI_IMPL_INST().at(index); }

  // The number of elements in the field.
  unsigned numElements() const { return VVI_IMPL_INST().numElements(); }
  unsigned numInternalElements() const { return VVI_IMPL_INST().numInternalElements(); }
  unsigned numGhostElements() const { return VVI_IMPL_INST().numGhostElements(); }
  unsigned size() const { return VVI_IMPL_INST().size(); }

  // Zero out the field elements.
  void Zero() { return VVI_IMPL_INST().Zero(); }

  // Methods to apply limits to Field data members.
  void applyMin(const DataType& dataMin) { VVI_IMPL_INST().applyMin(dataMin); }
  void applyMax(const DataType& dataMax) { VVI_IMPL_INST().applyMax(dataMax); }

  void applyScalarMin(const double dataMin) { VVI_IMPL_INST().applyScalarMin(dataMin); }
  void applyScalarMax(const double dataMax) { VVI_IMPL_INST().applyScalarMax(dataMax); }

  // Standard field additive operators.
  Field operator+(const Field& rhs) const { return VVI_IMPL_INST().operator+(rhs); }
  Field operator-(const Field& rhs) const { return VVI_IMPL_INST().operator-(rhs); }

  Field& operator+=(const Field& rhs) { return VVI_IMPL_INST().operator+=(rhs); }
  Field& operator-=(const Field& rhs) { return VVI_IMPL_INST().operator-=(rhs); }

  Field operator+(const DataType& rhs) const { return VVI_IMPL_INST().operator+(rhs); }
  Field operator-(const DataType& rhs) const { return VVI_IMPL_INST().operator-(rhs); }

  Field& operator+=(const DataType& rhs) { return VVI_IMPL_INST().operator+=(rhs); }
  Field& operator-=(const DataType& rhs) { return VVI_IMPL_INST().operator-=(rhs); }


  // Multiplication and division by scalar(s)
  Field operator*(const Field& rhs) const { return VVI_IMPL_INST().operator*(rhs); }
  Field operator/(const Field& rhs) const { return VVI_IMPL_INST().operator/(rhs); }

  Field& operator*=(const Field& rhs) { VVI_IMPL_INST().operator*=(rhs); return *this; }
  Field& operator/=(const Field& rhs) { VVI_IMPL_INST().operator/=(rhs); return *this; }

  Field operator*(const Scalar& rhs) const { return VVI_IMPL_INST().operator*(rhs); }
  Field operator/(const Scalar& rhs) const { return VVI_IMPL_INST().operator/(rhs); }

  Field& operator*=(const Scalar& rhs) { VVI_IMPL_INST().operator*=(rhs); return *this; }
  Field& operator/=(const Scalar& rhs) { VVI_IMPL_INST().operator/=(rhs); return *this; }
  
  // Some useful reduction operations.
  DataType sumElements() const { return VVI_IMPL_INST().sumElements(); }
  DataType min() const { return VVI_IMPL_INST().min(); }
  DataType max() const { return VVI_IMPL_INST().max(); }

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const { return VVI_IMPL_INST().localSumElements(); }
  DataType localMin() const { return VVI_IMPL_INST().localMin(); }
  DataType localMax() const { return VVI_IMPL_INST().localMax(); }

  // Comparison operators (Field-Field element wise).
  bool operator==(const Field& rhs) const { return VVI_IMPL_INST().operator==(rhs); }
  bool operator!=(const Field& rhs) const { return VVI_IMPL_INST().operator!=(rhs); }
  bool operator>(const Field& rhs) const { return VVI_IMPL_INST().operator>(rhs); }
  bool operator<(const Field& rhs) const { return VVI_IMPL_INST().operator<(rhs); }
  bool operator>=(const Field& rhs) const { return VVI_IMPL_INST().operator>=(rhs); }
  bool operator<=(const Field& rhs) const { return VVI_IMPL_INST().operator<=(rhs); }

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const { return VVI_IMPL_INST().operator==(rhs); }
  bool operator!=(const DataType& rhs) const { return VVI_IMPL_INST().operator!=(rhs); }
  bool operator>(const DataType& rhs) const { return VVI_IMPL_INST().operator>(rhs); }
  bool operator<(const DataType& rhs) const { return VVI_IMPL_INST().operator<(rhs); }
  bool operator>=(const DataType& rhs) const { return VVI_IMPL_INST().operator>=(rhs); }
  bool operator<=(const DataType& rhs) const { return VVI_IMPL_INST().operator<=(rhs); }

  // Test if this Field is in a valid, internally consistent state.
  bool valid() const { return VVI_IMPL_INST().valid(); }

  // Provide the standard iterator methods over the field.
  iterator begin() { return VVI_IMPL_INST().begin(); }
  iterator end() { return VVI_IMPL_INST().end(); }
  iterator internalBegin() { return VVI_IMPL_INST().internalBegin(); }
  iterator internalEnd() { return VVI_IMPL_INST().internalEnd(); }
  iterator ghostBegin() { return VVI_IMPL_INST().ghostBegin(); }
  iterator ghostEnd() { return VVI_IMPL_INST().ghostEnd(); }

  const_iterator begin() const { return VVI_IMPL_INST().begin(); }
  const_iterator end() const { return VVI_IMPL_INST().end(); }
  const_iterator internalBegin() const { return VVI_IMPL_INST().internalBegin(); }
  const_iterator internalEnd() const { return VVI_IMPL_INST().internalEnd(); }
  const_iterator ghostBegin() const { return VVI_IMPL_INST().ghostBegin(); }
  const_iterator ghostEnd() const { return VVI_IMPL_INST().ghostEnd(); }

  // Index operator.
  DataType& operator[](const unsigned int index) { return VVI_IMPL_INST().operator[](index); }
  const DataType& operator[](const unsigned int index) const { return VVI_IMPL_INST().operator[](index); }

  // Required functions from FieldBase
  void setNodeList(const NodeList<Dimension>& nodeList) { VVI_IMPL_INST().setNodeList(nodeList); }
  void resizeField(unsigned size) { VVI_IMPL_INST().resizeField(size); }
  void resizeFieldInternal(unsigned size, unsigned oldFirstGhostNode) { VVI_IMPL_INST().resizeFieldInternal(size, oldFirstGhostNode); }
  void resizeFieldGhost(unsigned size) { VVI_IMPL_INST().resizeFieldGhost(size); }
  void deleteElement(int nodeID) { VVI_IMPL_INST().deleteElement(nodeID); }
  void deleteElements(const std::vector<int>& nodeIDs) { VVI_IMPL_INST().deleteElements(nodeIDs); }

  std::vector<char> packValues(const std::vector<int>& nodeIDs) const { return VVI_IMPL_INST().packValues(nodeIDs); }
  void unpackValues(const std::vector<int>& nodeIDs, const std::vector<char>& buffer) { VVI_IMPL_INST().unpackValues(nodeIDs, buffer); }

  void copyElements(const std::vector<int>& fromIndices, const std::vector<int>& toIndices) { VVI_IMPL_INST().copyElements(fromIndices, toIndices); }

  bool fixedSizeDataType() const { return VVI_IMPL_INST().fixedSizeDataType(); }
  int numValsInDataType() const { return VVI_IMPL_INST().numValsInDataType(); }
  int sizeofDataType() const { return VVI_IMPL_INST().sizeofDataType(); }

  int computeCommBufferSize(const std::vector<int>& packIndices, const int sendProc, const int recvProc) const 
    { return VVI_IMPL_INST().computeCommBufferSize(packIndices, sendProc, recvProc); }

  // Serialization methods
  std::vector<char> serialize() const { return VVI_IMPL_INST().serialize(); }
  void deserialize(const std::vector<char>& buf) { VVI_IMPL_INST().deserialize(buf); }

  // Provide std::vector copies of the data.  This is mostly useful for the
  // python interface.
  std::vector<DataType> internalValues() const { return VVI_IMPL_INST().internalValues(); }
  std::vector<DataType> ghostValues() const { return VVI_IMPL_INST().ghostValues(); }
  std::vector<DataType> allValues() const { return VVI_IMPL_INST().allValues(); }

  // Functions to help with storing the field in a Sidre datastore.
  axom::sidre::DataTypeId getAxomTypeID() const { return VVI_IMPL_INST().getAxomTypeID(); }

  //FieldBase Forwarding Functions...
  const NodeList<Dimension>& nodeList() const { return VVI_IMPL_INST().nodeList(); }
  const NodeList<Dimension>* nodeListPtr() const { return VVI_IMPL_INST().nodeListPtr(); }
  ));
  
#endif // !defined(SPHERAL_ENABLE_VVI)

} // namespace Spheral

#include "FieldInline.hh"

#endif
