//---------------------------------Spheral++----------------------------------//
// FieldSpan -- provides a reference view (span) of the elements in an existing
// Field.
//
// Created by JMO, Mon Apr 28 15:05:15 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldSpan__
#define __Spheral_FieldSpan__

#include <span>

namespace Spheral {

template<typename Dimension, typename DataType> class Field;

template<typename Dimension, typename DataType>
class FieldSpan {
   
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;

  using FieldDimension = Dimension;
  using FieldDataType = DataType;
  using value_type = DataType;      // STL compatibility.

  using iterator = typename std::span<DataType>::iterator;
  // using const_iterator = typename std::span<DataType>::const_iterator;  // Not until C++23

  // Constructors, destructor
  FieldSpan(Field<Dimension, DataType>& field);
  FieldSpan(FieldSpan& rhs) = default;
  FieldSpan(FieldSpan&& rhs) = default;
  virtual ~FieldSpan() = default;

  // Assignment
  FieldSpan& operator=(FieldSpan& rhs) = default;
  FieldSpan& operator=(const DataType& rhs);

  // Element access.
  DataType& operator()(size_t index);
  const DataType& operator()(size_t index) const;

  DataType& at(size_t index);
  const DataType& at(size_t index) const;

  DataType& operator[](const size_t index);
  const DataType& operator[](const size_t index) const;

  // The number of elements in the field.
  size_t numElements()         const { return mDataSpan.size(); }
  size_t numInternalElements() const { return mNumInternalElements; }
  size_t numGhostElements()    const { return mNumGhostElements; }

  // Methods to apply limits to Field data members.
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const Scalar& dataMin);
  void applyScalarMax(const Scalar& dataMax);

  // Standard field additive operators.
  FieldSpan& operator+=(const FieldSpan& rhs);
  FieldSpan& operator-=(const FieldSpan& rhs);

  FieldSpan& operator+=(const DataType& rhs);
  FieldSpan& operator-=(const DataType& rhs);

  // Multiplication and division by scalar(s)
  FieldSpan& operator*=(const FieldSpan<Dimension, Scalar>& rhs);
  FieldSpan& operator/=(const FieldSpan<Dimension, Scalar>& rhs);

  FieldSpan& operator*=(const Scalar& rhs);
  FieldSpan& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const;
  DataType localMin() const;
  DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  bool operator==(const FieldSpan& rhs) const;
  bool operator!=(const FieldSpan& rhs) const;
  // bool operator> (const FieldSpan& rhs) const;
  // bool operator< (const FieldSpan& rhs) const;
  // bool operator>=(const FieldSpan& rhs) const;
  // bool operator<=(const FieldSpan& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator> (const DataType& rhs) const;
  bool operator< (const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  // Provide the standard iterator methods over the field.
  iterator begin() const                                              { return mDataSpan.begin(); }
  iterator end() const                                                { return mDataSpan.end(); }
  iterator internalBegin() const                                      { return mDataSpan.begin(); }
  iterator internalEnd() const                                        { return mDataSpan.begin() + mNumInternalElements; }
  iterator ghostBegin() const                                         { return mDataSpan.begin() + mNumInternalElements; }
  iterator ghostEnd() const                                           { return mDataSpan.end(); }

  // No default constructor.
  FieldSpan() = delete;

protected:
  //--------------------------- Protected Interface ---------------------------//
  // Private Data
  std::span<DataType> mDataSpan;
  size_t mNumInternalElements, mNumGhostElements;
};

} // namespace Spheral

#include "FieldSpanInline.hh"

#endif
