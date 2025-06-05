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

#include "Field/FieldSpanListBase.hh"

#include <span>

namespace Spheral {

// Forward declarations.
template<typename Dimension, typename DataType> class FieldSpan;

template<typename Dimension, typename DataType>
class FieldSpanList: public FieldSpanListBase<Dimension> {
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
  FieldSpanList(std::span<FieldSpan<Dimension, DataType>>& rhs);
  FieldSpanList(FieldSpanList& rhs) = default;
  FieldSpanList(FieldSpanList&& rhs) = default;
  virtual ~FieldSpanList() = default;

  // Assignment
  FieldSpanList& operator=(FieldSpanList& rhs) = default;
  FieldSpanList& operator=(const DataType& rhs);

  // Force the FieldSpan members of this FieldSpanList to be equal to those of
  // another FieldSpanList.
  void assignFields(FieldSpanList& rhs);

  // Provide the standard iterators over the FieldSpans
  iterator begin()                 { return mFieldSpans.begin(); }
  iterator end()                   { return mFieldSpans.end(); }
  reverse_iterator rbegin()        { return mFieldSpans.rbegin(); }
  reverse_iterator rend()          { return mFieldSpans.rend(); }

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

  FieldSpanList<Dimension, DataType>& operator*=(const FieldSpanList<Dimension, Scalar>& rhs);
  FieldSpanList<Dimension, DataType>& operator*=(const Scalar& rhs);

  FieldSpanList<Dimension, DataType>& operator/=(const FieldSpanList<Dimension, Scalar>& rhs);
  FieldSpanList<Dimension, DataType>& operator/=(const Scalar& rhs);

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
  size_t numFields() const      { return mFieldSpans.size(); }
  size_t size() const           { return mFieldSpans.size(); }
  bool empty() const            { return mFieldSpans.empty(); }

  // The number of nodes in the FieldSpanList.
  size_t numElements() const;
  
  // The number of internal nodes in the FieldSpanList.
  size_t numInternalElements() const;
  
  // The number of ghost nodes in the FieldSpanList.
  size_t numGhostElements() const;

  // No default constructor.
  FieldSpanList() = delete;

private:
  //--------------------------- Private Interface ---------------------------//
  StorageType mFieldSpans;
};

}

#include "FieldSpanListInline.hh"

#endif
