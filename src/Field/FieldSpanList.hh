//---------------------------------Spheral++----------------------------------//
// FieldSpanList -- A list container for FieldSpans
//
// Mostly a thin reference container for FieldSpans corresponding to the Fields
// in a normal FieldList.
//
// Created by JMO, Thu May  1 15:20:11 PDT 2025
//----------------------------------------------------------------------------//
#ifndef __Spheral__FieldSpanList__
#define __Spheral__FieldSpanList__

#include <span>

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class FieldSpan;

template<typename Dimension, typename DataType>
class FieldSpanList {
public:
  //--------------------------- Public Interface ---------------------------//
  using Scalar = typename Dimension::Scalar;
  
  using FieldDimension = Dimension;
  using FieldDataType = DataType;

  using value_type = FieldSpan<Dimension, DataType>;    // STL compatibility
  using StorageType = std::span<value_type>;

  using iterator = typename StorageType::iterator;
  using reverse_iterator = typename StorageType::reverse_iterator;

  // Constructors, destructor
  FieldSpanList() = default;
  FieldSpanList(std::span<FieldSpan<Dimension, DataType>>& rhs);
  FieldSpanList(FieldSpanList& rhs) = default;
  FieldSpanList(FieldSpanList&& rhs) = default;
  virtual ~FieldSpanList() = default;

  // Assignment
  FieldSpanList& operator=(FieldSpanList& rhs) = default;
  FieldSpanList& operator=(const DataType& rhs);

  // Provide the standard iterators over the FieldSpans
  iterator begin()                 { return mSpanFieldSpans.begin(); }
  iterator end()                   { return mSpanFieldSpans.end(); }
  reverse_iterator rbegin()        { return mSpanFieldSpans.rbegin(); }
  reverse_iterator rend()          { return mSpanFieldSpans.rend(); }

  // Index operator.
  value_type& operator[](const size_t index);
  const value_type& operator[](const size_t index) const;

  value_type& at(const size_t index);
  const value_type& at(const size_t index) const;

  // Provide direct access to Field elements
  DataType& operator()(const size_t fieldIndex,
                       const size_t nodeIndex);
  const DataType& operator()(const size_t fieldIndex,
                             const size_t nodeIndex) const;

  // Reproduce the standard Field operators for FieldSpanLists.
  void applyMin(const DataType& dataMin);
  void applyMax(const DataType& dataMax);

  void applyScalarMin(const Scalar dataMin);
  void applyScalarMax(const Scalar dataMax);

  FieldSpanList& operator+=(const FieldSpanList& rhs);
  FieldSpanList& operator-=(const FieldSpanList& rhs);

  FieldSpanList& operator+=(const DataType& rhs);
  FieldSpanList& operator-=(const DataType& rhs);

  FieldSpanList& operator*=(const FieldSpanList<Dimension, Scalar>& rhs);
  FieldSpanList& operator*=(const Scalar& rhs);

  FieldSpanList& operator/=(const FieldSpanList<Dimension, Scalar>& rhs);
  FieldSpanList& operator/=(const Scalar& rhs);

  // Some useful reduction operations.
  DataType sumElements() const;
  DataType min() const;
  DataType max() const;

  // Some useful reduction operations (local versions -- no MPI reductions)
  DataType localSumElements() const;
  DataType localMin() const;
  DataType localMax() const;

  // Comparison operators (Field-Field element wise).
  bool operator==(const FieldSpanList& rhs) const;
  bool operator!=(const FieldSpanList& rhs) const;
  bool operator>(const FieldSpanList& rhs) const;
  bool operator<(const FieldSpanList& rhs) const;
  bool operator>=(const FieldSpanList& rhs) const;
  bool operator<=(const FieldSpanList& rhs) const;

  // Comparison operators (Field-value element wise).
  bool operator==(const DataType& rhs) const;
  bool operator!=(const DataType& rhs) const;
  bool operator>(const DataType& rhs) const;
  bool operator<(const DataType& rhs) const;
  bool operator>=(const DataType& rhs) const;
  bool operator<=(const DataType& rhs) const;

  // The number of fields in the FieldSpanList.
  size_t numFields() const      { return mSpanFieldSpans.size(); }
  size_t size() const           { return mSpanFieldSpans.size(); }
  bool empty() const            { return mSpanFieldSpans.empty(); }

  // The number of nodes in the FieldSpanList.
  size_t numElements() const;
  
  // The number of internal nodes in the FieldSpanList.
  size_t numInternalElements() const;
  
  // The number of ghost nodes in the FieldSpanList.
  size_t numGhostElements() const;

protected:
  //--------------------------- Protected Interface ---------------------------//
  StorageType mSpanFieldSpans;
};

}

#include "FieldSpanListInline.hh"

#endif
