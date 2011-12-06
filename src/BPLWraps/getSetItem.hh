#ifndef __Spheral_getSetItem_hh__
#define __Spheral_getSetItem_hh__

#include <iostream>

#ifndef __GCCXML__
#include "boost/python.hpp"

#include "DBC.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Handle python indexing, including negative numbers.
//------------------------------------------------------------------------------
inline
int
handleIndex(int index, int begin, int end) {
  if (index < 0) index += end;
  index = std::min(end, std::max(begin, index));
  if (index < begin || index > end) {
    std::cerr << "Indexing error: index out of range: "
	      << index << " [" << begin << ", " << end << "]"
	      << std::endl;
    throw("Indexing error");
  }
  return index;
}

//------------------------------------------------------------------------------
// Generic methods to handle __getitem__ and __setitem__ methods.
//------------------------------------------------------------------------------
template<typename ContainerType, typename ElementType>
ElementType&
getItemByReference(ContainerType& self,
                   int index) {
  index = handleIndex(index, 0, self.size());
  return self[index];
}

template<typename ContainerType, typename ElementType>
ElementType
getItemByValue(ContainerType& self,
               int index) {
  index = handleIndex(index, 0, self.size());
  return self[index];
}

template<typename ContainerType>
boost::python::list
getSlice(ContainerType& self,
	 int start,
	 int stop) {
  boost::python::list result;
  start = handleIndex(start, 0, self.size());
  stop = handleIndex(stop,  0, self.size());
  CHECK(start >= 0 && start < self.size());
  CHECK(stop >= start && stop <= self.size());
  for (typename ContainerType::iterator itr = self.begin() + start;
       itr < self.begin() + stop;
       ++itr) {
    result.append(*itr);
  }
  return result;
}

template<typename ContainerType, typename ElementType>
void
setItemByReference(ContainerType& self,
                   int index,
                   ElementType& value) {
  index = handleIndex(index, 0, self.size());
  self[index] = value;
}

template<typename ContainerType, typename ElementType>
void
setItemByValue(ContainerType& self,
               int index,
               ElementType value) {
  index = handleIndex(index, 0, self.size());
  self[index] = value;
}

template<typename ContainerType, typename ElementType>
void
setSlice(ContainerType& self,
	 int start,
	 int stop,
	 const boost::python::list& values) {
  start = handleIndex(start, 0, self.size());
  stop = handleIndex(stop,  0, self.size());
  CHECK(start >= 0 && start < self.size());
  CHECK(stop >= start && stop <= self.size());
  if (stop - start != values.attr("__len__")()) {
    std::cerr << "Error: incorrect size in __setslice__ operation on " << &self << std::endl;
  } else {
    for (int i = start; i < stop; ++i) {
      self[i] = boost::python::extract<ElementType>(values[i]);
    }
  }
}

}

#else

#include "fakeboost.hh"

// Forward declarations.
namespace Spheral {

  int handleIndex(int index, int begin, int end);

  template<typename ContainerType, typename ElementType>
  ElementType&
  getItemByReference(ContainerType& self,
                     int index);

  template<typename ContainerType, typename ElementType>
  ElementType
  getItemByValue(ContainerType& self,
                 int index);

  template<typename ContainerType>
  boost::python::list
  getSlice(ContainerType& self,
           int start,
           int stop);

  template<typename ContainerType, typename ElementType>
  void
  setItemByReference(ContainerType& self,
                     int index,
                     ElementType& value);

  template<typename ContainerType, typename ElementType>
  void
  setItemByValue(ContainerType& self,
                 int index,
                 ElementType value);

  template<typename ContainerType, typename ElementType>
  void
  setSlice(ContainerType& self,
           int start,
           int stop,
           const boost::python::list& values);
}

#endif

#endif

