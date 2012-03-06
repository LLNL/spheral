#ifndef __Spheral_STLVectorWrap__
#define __Spheral_STLVectorWrap__

#include <vector>
#include <sstream>
#include <string>

#include "DBC.hh"

// BPL includes
#include "boost/python.hpp"

namespace Spheral {

using namespace std;
using namespace boost::python;

//------------------------------------------------------------------------------
// Methods to access element wise into vectors.
//------------------------------------------------------------------------------
template<typename DataType>
inline
const DataType
getItemValue(const std::vector<DataType>& self, int index) {
  if (index < 0) index += self.size();
  REQUIRE(index < self.size());
  return self[index];
}

template<typename DataType>
inline
DataType&
getItemReference(std::vector<DataType>& self, int index) {
  if (index < 0) index += self.size();
  REQUIRE(index < self.size());
  return self[index];
}

template<typename DataType>
inline
void
setItem(std::vector<DataType>& self,
        int index, const DataType& value) {
  if (index < 0) index += self.size();
  REQUIRE(index < self.size());
  self[index] = value;
}

inline
int indexInRange(const int index, const int minValue, const int maxValue) {
  REQUIRE(index >= 0);
  REQUIRE(minValue >= 0);
  REQUIRE(maxValue >= 0);
  return max(minValue, min(maxValue, index));
}

//------------------------------------------------------------------------------
// Methods to support slice operations.
//------------------------------------------------------------------------------
template<typename ContainerType>
inline
ContainerType
getSlice(const ContainerType& self, 
         int low, int high) {
  REQUIRE(low >= 0);
  REQUIRE(high >= 0);
  low = indexInRange(low, 0, self.size());
  high = indexInRange(high, 0, self.size());
  const int sliceSize = max(0, high - low);
  ContainerType result(sliceSize);
  CHECK(result.size() == sliceSize);
  for (int index = low; index < high; ++index)
    result[index - low] = self[index];
  return result;
}

template<typename ContainerType>
inline
void
setSlice(ContainerType& self, 
         int low, int high, const ContainerType& sliceValues) {
  REQUIRE(low >= 0);
  REQUIRE(high >= 0);
  low = indexInRange(low, 0, self.size());
  high = indexInRange(high, 0, self.size());
  const int sliceSize = max(0, high - low);
  CHECK(sliceSize == sliceValues.size());
  for (int index = low; index < high; ++index)
    self[index] = sliceValues[index - low];
}

//------------------------------------------------------------------------------
// Specialized methods to access element wise into vectors of pointers.
//------------------------------------------------------------------------------
template<typename DataType>
inline
DataType*
getItemPtr(std::vector<DataType*>& self, const unsigned int index) {
  REQUIRE(index < self.size());
  return self[index];
}

template<typename DataType>
void
setItemPtr(std::vector<DataType*>& self,
           const unsigned int index, DataType& value) {
  REQUIRE(index < self.size());
  self[index] = &value;
}

//------------------------------------------------------------------------------
// Helper to resize a vector.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
resizeSTLVector(std::vector<DataType>& self, const unsigned int size) {
  self.resize(size);
}

//------------------------------------------------------------------------------
// Helper method to provide a nice string representation of a vector.
//------------------------------------------------------------------------------
template<typename DataType>
inline
std::string
printSTLVector(const std::vector<DataType>& self) {
  std::stringstream resultStream;
  resultStream << "vector(";
  for (typename std::vector<DataType>::const_iterator elementItr = self.begin();
       elementItr < self.end();
       ++elementItr) {
    resultStream << *elementItr;
    if (elementItr < self.end() - 1) resultStream << ", ";
  }
  resultStream << ")" << '\0';
  return resultStream.str();
}

template<typename DataType>
inline
std::string
printSTLVectorPtr(const std::vector<DataType*>& self) {
  std::stringstream resultStream;
  resultStream << "vector(";
  for (typename std::vector<DataType*>::const_iterator elementItr = self.begin();
       elementItr < self.end();
       ++elementItr) {
    resultStream << *elementItr;
    if (elementItr < self.end() - 1) resultStream << ", ";
  }
  resultStream << ")" << '\0';
  return resultStream.str();
}

//------------------------------------------------------------------------------
// Print a nice string representation for iterators.
//------------------------------------------------------------------------------
template<typename IteratorType>
inline
std::string
printSTLVectorIterator(const IteratorType& self) {
  std::stringstream resultStream;
  resultStream << "Iterator(" << &self << ")" << ends;
  return resultStream.str();
}

//------------------------------------------------------------------------------
// A templated method to provide a "standard" wrapper for STL vectors.
// This version returns elements by value.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
wrapSTLVectorByValue(const std::string label) {

  typedef std::vector<DataType> VectorType;

  class_<VectorType>(label.c_str(), init<>())

    // Constructors.
    .def(init<unsigned int>())

    // Standard vector methods.
    .def("size", &VectorType::size)
    .def("resize", resizeSTLVector<DataType>)

    .def("__len__", &VectorType::size)
    .def("__getitem__", getItemValue<DataType>)
    .def("__setitem__", setItem<DataType>)

    .def("__getslice__", getSlice< std::vector<DataType> >)
    .def("__setslice__", setSlice< std::vector<DataType> >)

    // Iterator support.
    .def("__iter__", boost::python::iterator<VectorType>())

    // String representations.
    .def("__str__", printSTLVector<DataType>)
    .def("__repr__", printSTLVector<DataType>)

    ;
}

//------------------------------------------------------------------------------
// A templated method to provide a "standard" wrapper for STL vectors.
// This version returns elements by reference.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
wrapSTLVectorByReference(const std::string label) {

  typedef std::vector<DataType> VectorType;

  class_<VectorType>(label.c_str(), init<>())

    // Constructors.
    .def(init<unsigned int>())

    // Standard vector methods.
    .def("size", &VectorType::size)
    .def("resize", resizeSTLVector<DataType>)

    .def("__len__", &VectorType::size)
    .def("__getitem__", getItemReference<DataType>, return_internal_reference<>())
    .def("__setitem__", setItem<DataType>)

    .def("__getslice__", getSlice< std::vector<DataType> >)
    .def("__setslice__", setSlice< std::vector<DataType> >)

    // Iterator support.
    .def("__iter__", boost::python::iterator<VectorType, return_internal_reference<> >())

    // String representations.
    .def("__str__", printSTLVector<DataType>)
    .def("__repr__", printSTLVector<DataType>)

    ;
}

//------------------------------------------------------------------------------
// A templated method to provide a "standard" wrapper for STL vectors of 
// pointers.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
wrapSTLVectorPtr(const std::string label) {

  typedef std::vector<DataType*> VectorType;

  class_<VectorType>(label.c_str(), init<>())

    // Constructors.
    .def(init<unsigned int>())

    // Standard vector methods.
    .def("size", &VectorType::size)
    .def("resize", resizeSTLVector<DataType*>)

    .def("__len__", &VectorType::size)
    .def("__getitem__", getItemPtr<DataType>, return_internal_reference<>())
    .def("__setitem__", setItemPtr<DataType>)

    // Iterator support.
    .def("__iter__", boost::python::iterator<VectorType, return_internal_reference<> >())

    // String representations.
    .def("__str__", printSTLVectorPtr<DataType>)
    .def("__repr__", printSTLVectorPtr<DataType>)

    ;
}

//------------------------------------------------------------------------------
// Provide a standard wrapper for STL vector<>::iterators with pointers.
//------------------------------------------------------------------------------
template<typename DataType>
inline
void
wrapSTLVectorPtrIterator(const std::string label) {

  typedef std::vector<DataType*> VectorType;
  typedef typename std::vector<DataType*>::iterator IteratorType;
  typedef typename std::vector<DataType*>::const_iterator ConstIteratorType;
  
  // Export the iterator.
  class_<IteratorType>(label.c_str())
    // String representations.
    .def("__str__", printSTLVectorIterator<IteratorType>)
    .def("__repr__", printSTLVectorPtr<IteratorType>)
    ;

  // Export the const_iterator.
  class_<ConstIteratorType>(("const_" + label).c_str())
    // String representations.
    .def("__str__", printSTLVectorIterator<ConstIteratorType>)
    .def("__repr__", printSTLVectorPtr<ConstIteratorType>)
    ;

}

}

#endif
