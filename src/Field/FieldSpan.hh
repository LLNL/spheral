//---------------------------------Spheral++----------------------------------//
// FieldSpan -- provides a reference view (span) of the elements in an existing
// Field.
//
// Created by JMO, Mon Apr 28 15:05:15 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral_FieldSpan__
#define __Spheral_FieldSpan__

#include "FieldSpanBase.hh"

#include <span>

namespace Spheral {

template<typename Dimension, typename DataType>
class FieldSpan: public FieldSpanBase<Dimension> {
   
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;

  using FieldDimension = Dimension;
  using FieldDataType = DataType;
  using value = DataType;      // STL compatibility.

  using iterator = typename std::span<DataType>::iterator;
  // using const_iterator = typename std::span<DataType>::const_iterator;  // Not until C++23

  // Constructors, destructor
  FieldSpan(Field<Dimension, DataType>& field);
  virtual ~FieldSpan() = default;

  // Assignment operator.
  // virtual FieldSpanBase<Dimension>& operator=(FieldSpanBase<Dimension>& rhs) override;
  FieldSpan(FieldSpan& rhs) = default;
  FieldSpan& operator=(const DataType& rhs);

  // Required method to test equivalence with a FieldSpanBase.
  virtual bool operator==(const FieldSpanBase<Dimension>& rhs) const override;

  // Element access.
  DataType& operator()(size_t index);
  const DataType& operator()(size_t index) const;

  DataType& at(size_t index);
  const DataType& at(size_t index) const;

  DataType& operator[](const size_t index);
  const DataType& operator[](const size_t index) const;

  // The number of elements in the field.
  size_t numElements() const;
  size_t numInternalElements() const;
  size_t numGhostElements() const;
  virtual size_t size() const override;

  // Zero out the field elements.
  virtual void Zero() override;

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
  bool operator> (const FieldSpan& rhs) const;
  bool operator< (const FieldSpan& rhs) const;
  bool operator>=(const FieldSpan& rhs) const;
  bool operator<=(const FieldSpan& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator> (const DataType& rhs) const;
  bool operator< (const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  // Provide the standard iterator methods over the field.
  iterator begin();
  iterator end();
  iterator internalBegin();
  iterator internalEnd();
  iterator ghostBegin();
  iterator ghostEnd();

  // No default constructor.
  FieldSpan() = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  // Private Data
  std::span<DataType> mDataSpan;
  size_t mNumInternalElements, mNumGhostElements;
};

} // namespace Spheral

#include "FieldSpanInline.hh"

#endif
