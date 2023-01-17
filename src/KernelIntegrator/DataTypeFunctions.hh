//---------------------------------Spheral++----------------------------------//
// DataTypeFunctions
//
// Allows for math operations on templated data. The size data is to allow for
// size checking, and doesn't necessarily represent the total number of primitives. 
//----------------------------------------------------------------------------//
#ifndef __Spheral_DataTypeFunctions_hh__
#define __Spheral_DataTypeFunctions_hh__

// #include <algorithm> // std::transform
#include <tuple>
#include <utility> // std::pair
#include <vector>

#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Utilities/DataTypeTraits.hh"

namespace Spheral {

// General definition: applies to most types defined in DataTypeTraits
template<typename DataType>
struct DataTypeFunctions {
  static int size(const DataType& x) { return DataTypeTraits<DataType>::numElements(x); }
  static void initialize(DataType& x) { /* do nothing */ }
  static void zeroOutData(DataType& x) { x = DataTypeTraits<DataType>::zero(); }
  static void addToData(DataType& x, const DataType& y) { x += y; }
  static void subtractFromData(DataType& x, const DataType& y) { x -= y; }
  static void multiplyData(DataType& x, const double y) { x *= y; }
  static void copyData(DataType& x, const DataType& y) { x = y; }
}; 

// vector<ElementType>
template<typename ElementType>
struct DataTypeFunctions<std::vector<ElementType>> {
  typedef std::vector<ElementType> DataType;
  static void initialize(DataType& x) { /* do nothing */ }
  static int size(const DataType& x) { return x.size(); }
  static void zeroOutData(DataType& x) {
    std::fill(x.begin(), x.end(), DataTypeTraits<ElementType>::zero());
  }
  static void addToData(DataType& x, const DataType& y) {
    auto size = x.size();
    CHECK(size == y.size());
    for (auto i = 0u; i < size; ++i) {
      DataTypeFunctions<ElementType>::addToData(x[i], y[i]);
    }
    // std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::plus<ElementType>());
  }
  static void subtractFromData(DataType& x, const DataType& y) {
    auto size = x.size();
    CHECK(size == y.size());
    for (auto i = 0u; i < size; ++i) {
      DataTypeFunctions<ElementType>::addToData(x[i], y[i]);
    }
    // std::transform(x.begin(), x.end(), y.begin(), x.begin(), std::minus<ElementType>());
  }
  static void multiplyData(DataType& x, const double y) {
    auto size = x.size();
    for (auto i = 0u; i < size; ++i) {
      DataTypeFunctions<ElementType>::multiplyData(x[i], y);
    }
  }
  static void copyData(DataType& x, const DataType& y) { x = y; }
};

// Field<Dimension, ElementType>
template<typename Dimension, typename ElementType>
struct DataTypeFunctions<Field<Dimension, ElementType>> {
  typedef Field<Dimension, ElementType> DataType;
  static void initialize(DataType& x) { /* do nothing */ }
  static int size(const DataType& x) { return x.size(); }
  static void zeroOutData(DataType& x) {
    // x = 0;
    for (auto it = x.begin(); it != x.end(); ++it) {
      DataTypeFunctions<ElementType>::zeroOutData(*it);
    }
  }
  static void addToData(DataType& x, const DataType& y) {
    CHECK(x.size() == y.size());
    // x += y;
    for (auto it = std::make_pair(x.begin(), y.begin());
         it.first != x.end(); ++it.first, ++it.second) {
      DataTypeFunctions<ElementType>::addToData(*it.first, *it.second);
    }
  }
  static void subtractFromData(DataType& x, const DataType& y) {
    CHECK(x.size() == y.size());
    // x -= y;
    for (auto it = std::make_pair(x.begin(), y.begin());
         it.first != x.end(); ++it.first, ++it.second) {
      DataTypeFunctions<ElementType>::subtractFromData(*it.first, *it.second);
    }
  }
  static void multiplyData(DataType& x, const double y) {
    for (auto it = x.begin(); it != x.end(); ++it) {
      DataTypeFunctions<ElementType>::multiplyData(*it, y);
    }
  }
  static void copyData(DataType& x, const DataType& y) { x = y; }
};
  
// FieldList<Dimension, ElementType>
template<typename Dimension, typename ElementType>
struct DataTypeFunctions<FieldList<Dimension, ElementType>> {
  typedef FieldList<Dimension, ElementType> DataType;
  typedef Field<Dimension, ElementType> FieldType;
  static void initialize(DataType& x) { x.copyFields(); }
  static int size(const DataType& x) { return x.size(); }
  static void zeroOutData(DataType& x) {
    // x = 0;
    for (auto it = x.begin(); it != x.end(); ++it) {
      DataTypeFunctions<FieldType>::zeroOutData(**it);
    }
  }
  static void addToData(DataType& x, const DataType& y) {
    CHECK(x.size() == y.size());
    // x += y;
    for (auto it = std::make_pair(x.begin(), y.begin());
         it.first != x.end(); ++it.first, ++it.second) {
      DataTypeFunctions<FieldType>::addToData(**it.first, **it.second);
    }
  }
  static void subtractFromData(DataType& x, const DataType& y) {
    CHECK(x.size() == y.size());
    // x -= y;
    for (auto it = std::make_pair(x.begin(), y.begin());
         it.first != x.end(); ++it.first, ++it.second) {
      DataTypeFunctions<FieldType>::subtractFromData(**it.first, **it.second);
    }
  }
  static void multiplyData(DataType& x, const double y) {
    for (auto it = x.begin(); it != x.end(); ++it) {
      DataTypeFunctions<FieldType>::multiplyData(**it, y);
    }
  }
  static void copyData(DataType& x, const DataType& y) {
    if (x.storageType() == FieldStorageType::ReferenceFields) {
      CHECK(x.size() == y.size());
      x.assignFields(y);
    }
    else {
      x = y;
    }
  }
};

} // end namespace Spheral

#endif
