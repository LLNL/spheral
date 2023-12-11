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
#ifndef __Spheral_FieldView_hh__
#define __Spheral_FieldView_hh__

#include "FieldBase.hh"
#include "axom/sidre.hpp"
#include "SphArray.hh"
#include "config.hh"

#include <vector>

#ifdef USE_UVM
#include "uvm_allocator.hh"
#endif

namespace Spheral {

template<typename Dimension, typename DataType>
class Field; 

template<typename Dimension, typename DataType>
class FieldView: 
    public chai::CHAICopyable{

public:
  //--------------------------- Public Interface ---------------------------//
  
  using FieldType = Field<Dimension, DataType>;

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::SymTensor SymTensor;
  
  typedef Dimension FieldDimension;
  typedef DataType FieldDataType;
  typedef DataType value_type;      // STL compatibility.

  using ContainerType = ManagedVector<DataType>;

  using iterator = typename ContainerType::iterator;
  using const_iterator = typename ContainerType::const_iterator;

  // Constructors.
  SPHERAL_HOST_DEVICE FieldView();
  SPHERAL_HOST FieldView(FieldType const& field);
  SPHERAL_HOST_DEVICE FieldView(const FieldView& field);

  // Destructor.
  SPHERAL_HOST_DEVICE ~FieldView();

  SPHERAL_HOST FieldType& operator*() const;
  SPHERAL_HOST FieldType* operator->() const;

  // Assignment operator.
  SPHERAL_HOST_DEVICE FieldView& operator=(const FieldView& rhs);
  SPHERAL_HOST_DEVICE FieldView& operator=(const ContainerType& rhs);
  SPHERAL_HOST FieldView& operator=(const DataType& rhs);

  // Element access.
  SPHERAL_HOST_DEVICE DataType& operator()(int index);
  SPHERAL_HOST_DEVICE DataType& operator()(int index) const;

  SPHERAL_HOST_DEVICE DataType& at(int index);
  SPHERAL_HOST_DEVICE const DataType& at(int index) const;

  // Index operator.
  SPHERAL_HOST_DEVICE DataType& operator[](const unsigned int index);
  SPHERAL_HOST_DEVICE DataType& operator[](const unsigned int index) const;

  // The number of elements in the field.
  SPHERAL_HOST_DEVICE unsigned numElements() const;
  SPHERAL_HOST_DEVICE virtual unsigned size() const;

  // Zero out the field elements.
  SPHERAL_HOST virtual void Zero();

  // Methods to apply limits to FieldView data members.
  SPHERAL_HOST void applyMin(const DataType& dataMin);
  SPHERAL_HOST void applyMax(const DataType& dataMax);

  SPHERAL_HOST void applyScalarMin(const double dataMin);
  SPHERAL_HOST void applyScalarMax(const double dataMax);

  // Standard field additive operators.
  SPHERAL_HOST FieldView operator+(const FieldView& rhs) const;
  SPHERAL_HOST FieldView operator-(const FieldView& rhs) const;

  SPHERAL_HOST FieldView& operator+=(const FieldView& rhs);
  SPHERAL_HOST FieldView& operator-=(const FieldView& rhs);

  SPHERAL_HOST FieldView operator+(const DataType& rhs) const;
  SPHERAL_HOST FieldView operator-(const DataType& rhs) const;

  SPHERAL_HOST FieldView& operator+=(const DataType& rhs);
  SPHERAL_HOST FieldView& operator-=(const DataType& rhs);

  // Multiplication and division by scalar(s)
  SPHERAL_HOST FieldView operator*(const FieldView<Dimension, Scalar>& rhs) const;
  SPHERAL_HOST FieldView operator/(const FieldView<Dimension, Scalar>& rhs) const;

  SPHERAL_HOST FieldView& operator*=(const FieldView<Dimension, Scalar>& rhs);
  SPHERAL_HOST FieldView& operator/=(const FieldView<Dimension, Scalar>& rhs);

  SPHERAL_HOST FieldView operator*(const Scalar& rhs) const;
  SPHERAL_HOST FieldView operator/(const Scalar& rhs) const;

  SPHERAL_HOST FieldView& operator*=(const Scalar& rhs);
  SPHERAL_HOST FieldView& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  SPHERAL_HOST DataType sumElements() const;
  SPHERAL_HOST DataType min() const;
  SPHERAL_HOST DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  SPHERAL_HOST DataType localSumElements() const;
  SPHERAL_HOST DataType localMin() const;
  SPHERAL_HOST DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  SPHERAL_HOST bool operator==(const FieldView& rhs) const;
  SPHERAL_HOST bool operator!=(const FieldView& rhs) const;
  SPHERAL_HOST bool operator>(const FieldView& rhs) const;
  SPHERAL_HOST bool operator<(const FieldView& rhs) const;
  SPHERAL_HOST bool operator>=(const FieldView& rhs) const;
  SPHERAL_HOST bool operator<=(const FieldView& rhs) const;

  // Comparison operators (Field-value element wise).
  SPHERAL_HOST bool operator==(const DataType& rhs) const;
  SPHERAL_HOST bool operator!=(const DataType& rhs) const;
  SPHERAL_HOST bool operator>(const DataType& rhs) const;
  SPHERAL_HOST bool operator<(const DataType& rhs) const;
  SPHERAL_HOST bool operator>=(const DataType& rhs) const;
  SPHERAL_HOST bool operator<=(const DataType& rhs) const;

  // Provide the standard iterator methods over the field.
  SPHERAL_HOST iterator begin();
  SPHERAL_HOST iterator end();
  SPHERAL_HOST const_iterator begin() const;
  SPHERAL_HOST const_iterator end() const;

  SPHERAL_HOST void move(chai::ExecutionSpace space, bool touch = true) const;

  // Required interface to make FieldView CHAICopyable.
  SPHERAL_HOST_DEVICE FieldView& operator=(std::nullptr_t) { mDataArray=nullptr; return *this; }
  SPHERAL_HOST_DEVICE void shallowCopy(const FieldView& field) {
    FieldView result(field);
    result.mDataArray.shallowCopy(field.mDataArray);
    *this = result;
  }

protected:
  //--------------------------- Protected Interface ---------------------------//
  ContainerType mDataArray;
  FieldType* mFieldPtr = nullptr;

};

} // namespace Spheral

#include "FieldViewInline.hh"

#else

// Forward declare the FieldView class.
namespace Spheral {
  template<typename Dimension, typename DataType> class FieldView;
} // namespace Spheral

#endif // __Spheral_FieldView_hh__
