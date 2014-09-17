#ifndef __PBGWRAPS_CXXTYPES__
#define __PBGWRAPS_CXXTYPES__

#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include "Geometry/Dimension.hh"
#include "Utilities/DataTypeTraits.hh"

using namespace std;

// namespace std {

//------------------------------------------------------------------------------
// Names!
//------------------------------------------------------------------------------
typedef pair<unsigned, unsigned> pair_unsigned_unsigned;
typedef pair<uint64_t, uint64_t> pair_ULL_ULL;
typedef pair<double, double> pair_double_double;
typedef pair<double, string> pair_double_string;
typedef pair<string, string> pair_string_string;

typedef pair<Spheral::Vector1d, Spheral::Vector1d> pair_Vector1d_Vector1d;
typedef pair<Spheral::Vector2d, Spheral::Vector2d> pair_Vector2d_Vector2d;
typedef pair<Spheral::Vector3d, Spheral::Vector3d> pair_Vector3d_Vector3d;

typedef pair<Spheral::Tensor1d, Spheral::Tensor1d> pair_Tensor1d_Tensor1d;
typedef pair<Spheral::Tensor2d, Spheral::Tensor2d> pair_Tensor2d_Tensor2d;
typedef pair<Spheral::Tensor3d, Spheral::Tensor3d> pair_Tensor3d_Tensor3d;

typedef vector<bool>               vector_of_bool;
typedef vector<char>               vector_of_char;
typedef vector<int>                vector_of_int;
typedef vector<float>              vector_of_float;
typedef vector<double>             vector_of_double;
typedef vector<string>             vector_of_string;
typedef vector<unsigned>           vector_of_unsigned;
typedef vector<uint64_t>           vector_of_ULL;
typedef vector<pair_unsigned_unsigned> vector_of_pair_unsigned_unsigned;
typedef vector<pair_ULL_ULL>       vector_of_pair_ULL_ULL;

typedef vector<bool>::iterator               vector_of_bool_iterator;
typedef vector<char>::iterator               vector_of_char_iterator;
typedef vector<int>::iterator                vector_of_int_iterator;
typedef vector<float>::iterator              vector_of_float_iterator;
typedef vector<double>::iterator             vector_of_double_iterator;
typedef vector<string>::iterator             vector_of_string_iterator;
typedef vector<unsigned>::iterator           vector_of_unsigned_iterator;
typedef vector<uint64_t>::iterator           vector_of_ULL_iterator;
typedef vector<pair_unsigned_unsigned>::iterator vector_of_pair_unsigned_unsigned_iterator;
typedef vector<pair_ULL_ULL>::iterator       vector_of_pair_ULL_ULL_iterator;

typedef vector<Spheral::Vector1d> vector_of_Vector1d;
typedef vector<Spheral::Vector2d> vector_of_Vector2d;
typedef vector<Spheral::Vector3d> vector_of_Vector3d;

typedef vector<Spheral::Tensor1d> vector_of_Tensor1d;
typedef vector<Spheral::Tensor2d> vector_of_Tensor2d;
typedef vector<Spheral::Tensor3d> vector_of_Tensor3d;

typedef vector<Spheral::SymTensor1d> vector_of_SymTensor1d;
typedef vector<Spheral::SymTensor2d> vector_of_SymTensor2d;
typedef vector<Spheral::SymTensor3d> vector_of_SymTensor3d;

typedef vector<Spheral::ThirdRankTensor1d> vector_of_ThirdRankTensor1d;
typedef vector<Spheral::ThirdRankTensor2d> vector_of_ThirdRankTensor2d;
typedef vector<Spheral::ThirdRankTensor3d> vector_of_ThirdRankTensor3d;

typedef vector<Spheral::Geom3Vector> vector_of_Geom3Vector;

typedef vector<Spheral::Box1d> vector_of_Box1d;
typedef vector<Spheral::GeomPolygon> vector_of_Polygon;
typedef vector<Spheral::GeomPolyhedron> vector_of_Polyhedron;

typedef vector<Spheral::Vector1d>::iterator vector_of_Vector1d_iterator;
typedef vector<Spheral::Vector2d>::iterator vector_of_Vector2d_iterator;
typedef vector<Spheral::Vector3d>::iterator vector_of_Vector3d_iterator;

typedef vector<Spheral::Tensor1d>::iterator vector_of_Tensor1d_iterator;
typedef vector<Spheral::Tensor2d>::iterator vector_of_Tensor2d_iterator;
typedef vector<Spheral::Tensor3d>::iterator vector_of_Tensor3d_iterator;

typedef vector<Spheral::SymTensor1d>::iterator vector_of_SymTensor1d_iterator;
typedef vector<Spheral::SymTensor2d>::iterator vector_of_SymTensor2d_iterator;
typedef vector<Spheral::SymTensor3d>::iterator vector_of_SymTensor3d_iterator;

typedef vector<Spheral::ThirdRankTensor1d>::iterator vector_of_ThirdRankTensor1d_iterator;
typedef vector<Spheral::ThirdRankTensor2d>::iterator vector_of_ThirdRankTensor2d_iterator;
typedef vector<Spheral::ThirdRankTensor3d>::iterator vector_of_ThirdRankTensor3d_iterator;

typedef vector<Spheral::Geom3Vector>::iterator vector_of_Geom3Vector_iterator;

typedef vector<Spheral::Box1d>::iterator vector_of_Box1d_iterator;
typedef vector<Spheral::GeomPolygon>::iterator vector_of_Polygon_iterator;
typedef vector<Spheral::GeomPolyhedron>::iterator vector_of_Polyhedron_iterator;

typedef vector<vector_of_char>        vector_of_vector_of_char;
typedef vector<vector_of_unsigned>    vector_of_vector_of_unsigned;
typedef vector<vector_of_int>         vector_of_vector_of_int;
typedef vector<vector_of_float>       vector_of_vector_of_float;
typedef vector<vector_of_double>      vector_of_vector_of_double;
typedef vector<vector_of_string>      vector_of_vector_of_string;

typedef vector<vector_of_Vector1d>       vector_of_vector_of_Vector1d;
typedef vector<vector_of_Vector2d>       vector_of_vector_of_Vector2d;
typedef vector<vector_of_Vector3d>       vector_of_vector_of_Vector3d;

typedef vector<vector_of_Tensor1d>       vector_of_vector_of_Tensor1d;
typedef vector<vector_of_Tensor2d>       vector_of_vector_of_Tensor2d;
typedef vector<vector_of_Tensor3d>       vector_of_vector_of_Tensor3d;

typedef vector<vector_of_SymTensor1d>       vector_of_vector_of_SymTensor1d;
typedef vector<vector_of_SymTensor2d>       vector_of_vector_of_SymTensor2d;
typedef vector<vector_of_SymTensor3d>       vector_of_vector_of_SymTensor3d;

typedef vector<vector_of_ThirdRankTensor1d>       vector_of_vector_of_ThirdRankTensor1d;
typedef vector<vector_of_ThirdRankTensor2d>       vector_of_vector_of_ThirdRankTensor2d;
typedef vector<vector_of_ThirdRankTensor3d>       vector_of_vector_of_ThirdRankTensor3d;

typedef vector<vector_of_vector_of_char>        vector_of_vector_of_vector_of_char;
typedef vector<vector_of_vector_of_unsigned>    vector_of_vector_of_vector_of_unsigned;
typedef vector<vector_of_vector_of_int>         vector_of_vector_of_vector_of_int;
typedef vector<vector_of_vector_of_float>       vector_of_vector_of_vector_of_float;
typedef vector<vector_of_vector_of_double>      vector_of_vector_of_vector_of_double;
typedef vector<vector_of_vector_of_string>      vector_of_vector_of_vector_of_string;

typedef vector<vector_of_char>::iterator        vector_of_vector_of_char_iterator;
typedef vector<vector_of_unsigned>::iterator    vector_of_vector_of_unsigned_iterator;
typedef vector<vector_of_int>::iterator         vector_of_vector_of_int_iterator;
typedef vector<vector_of_float>::iterator       vector_of_vector_of_float_iterator;
typedef vector<vector_of_double>::iterator      vector_of_vector_of_double_iterator;
typedef vector<vector_of_string>::iterator      vector_of_vector_of_string_iterator;

typedef vector<vector_of_Vector1d>::iterator       vector_of_vector_of_Vector1d_iterator;
typedef vector<vector_of_Vector2d>::iterator       vector_of_vector_of_Vector2d_iterator;
typedef vector<vector_of_Vector3d>::iterator       vector_of_vector_of_Vector3d_iterator;

typedef vector<vector_of_Tensor1d>::iterator       vector_of_vector_of_Tensor1d_iterator;
typedef vector<vector_of_Tensor2d>::iterator       vector_of_vector_of_Tensor2d_iterator;
typedef vector<vector_of_Tensor3d>::iterator       vector_of_vector_of_Tensor3d_iterator;

typedef vector<vector_of_SymTensor1d>::iterator       vector_of_vector_of_SymTensor1d_iterator;
typedef vector<vector_of_SymTensor2d>::iterator       vector_of_vector_of_SymTensor2d_iterator;
typedef vector<vector_of_SymTensor3d>::iterator       vector_of_vector_of_SymTensor3d_iterator;

typedef vector<vector_of_ThirdRankTensor1d>::iterator       vector_of_vector_of_ThirdRankTensor1d_iterator;
typedef vector<vector_of_ThirdRankTensor2d>::iterator       vector_of_vector_of_ThirdRankTensor2d_iterator;
typedef vector<vector_of_ThirdRankTensor3d>::iterator       vector_of_vector_of_ThirdRankTensor3d_iterator;

typedef vector<vector_of_vector_of_char>::iterator        vector_of_vector_of_vector_of_char_iterator;
typedef vector<vector_of_vector_of_unsigned>::iterator    vector_of_vector_of_vector_of_unsigned_iterator;
typedef vector<vector_of_vector_of_int>::iterator         vector_of_vector_of_vector_of_int_iterator;
typedef vector<vector_of_vector_of_float>::iterator       vector_of_vector_of_vector_of_float_iterator;
typedef vector<vector_of_vector_of_double>::iterator      vector_of_vector_of_vector_of_double_iterator;
typedef vector<vector_of_vector_of_string>::iterator      vector_of_vector_of_vector_of_string_iterator;

// }

namespace Spheral {

//------------------------------------------------------------------------------
// Index into a Container.
//------------------------------------------------------------------------------
template<typename Container>
inline
typename Container::value_type
indexContainer(Container& container,
               int index) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < n) {
    try {
      return container.at(index);
    } catch (std::out_of_range) {
      PyErr_SetString(PyExc_IndexError, "Container index out of range");
      return typename Container::value_type();
    }
  } else {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return typename Container::value_type();
  }
}

template<typename Container>
inline
typename Container::value_type*
indexContainerAsPointer(Container& container,
                        int index) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < n) {
    try {
      return &(container.at(index));
    } catch (std::out_of_range) {
      PyErr_SetString(PyExc_IndexError, "Container index out of range");
      return NULL;
    }
  } else {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return NULL;
  }
}

template<typename Container>
inline
typename Container::value_type&
indexContainerAsReference(Container& container,
                          int index) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < n) {
    try {
      return container.at(index);
    } catch (std::out_of_range) {
      PyErr_SetString(PyExc_IndexError, "Container index out of range");
      return typename Container::value_type();
    }
  } else {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return NULL;
  }
}

//------------------------------------------------------------------------------
// Extract slices from a container.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container
sliceContainer(Container& container,
               int index1,
               int index2) {
   const int n = DataTypeTraits<Container>::numElements(container);
   index1 = (index1 < 0 ? (n + 1 + index1) : index1);
   index2 = (index2 < 0 ? (n + 1 + index2) : index2);
   Container result;
   for (int i = index1; i < index2; ++i) result.push_back(indexContainer(container, i));
   return result;
}

//------------------------------------------------------------------------------
// Assign to a postion in a container.
//------------------------------------------------------------------------------
template<typename Container>
inline
int
assignToContainerIndex(Container& container, 
                       int index,
                       const typename Container::value_type& value) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < 0 or index >= container.size()) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return 1;
  } else {
    container[index] = value;
    return 0;
  }
}

template<typename Container>
inline
int
assignToContainerIndexPtr(Container& container, 
                          int index,
                          const typename Container::value_type& value) {
  const int n = DataTypeTraits<Container>::numElements(container);
  index = (index < 0 ? (n + 1 + index) : index);
  if (index < 0 or index >= container.size()) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return 1;
  } else {
    *(container[index]) = *(value);
    return 0;
  }
}

//------------------------------------------------------------------------------
// Assign to a slice in a container.
//------------------------------------------------------------------------------
template<typename Container>
inline
int
assignToSlice(Container& container, 
              int index1,
              int index2,
              const Container& values) {
   const int n = DataTypeTraits<Container>::numElements(container);
   index1 = (index1 < 0 ? (n + 1 + index1) : index1);
   index2 = (index2 < 0 ? (n + 1 + index2) : index2);
   const int nv = values.size();
   if (index2 - index1 != nv) {
     PyErr_SetString(PyExc_IndexError, "Container slices different sizes.");
     return -1;
   }
   for (int i = index1; i < index2; ++i) {
      if (assignToContainerIndex(container, i, values[i - index1]) != 0) return -1;
   }
   return 0;
}

//------------------------------------------------------------------------------
// Append a value to container of pointers.
//------------------------------------------------------------------------------
template<typename Container>
inline
int
appendToContainerOfPointers(Container& container, 
                            typename Container::value_type value) {
  container.push_back(value);
  return 0;
}

//------------------------------------------------------------------------------
// In-place concatenation.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container&
concatContainersInPlace(Container& lhs,
                        const Container& rhs) {
  const unsigned n = lhs.size() + rhs.size();
  std::copy(rhs.begin(), rhs.end(), std::back_inserter(lhs));
  assert(lhs.size() == n);
  return lhs;
}

//------------------------------------------------------------------------------
// In-place repeat.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container&
repeatContainerInPlace(Container& self,
                       const unsigned count) {
  const unsigned n0 = self.size();
  if (count == 0) {
    self = Container();
  } else {
    self.reserve(n0 * count);
    for (unsigned i = 0; i != count - 1; ++i) {
       std::copy(self.begin(), self.begin() + n0, std::back_inserter(self));
    }
  }
  assert(self.size() == n0 * count);
  return self;
}

//------------------------------------------------------------------------------
// Support concatenation of containers as addition.
//------------------------------------------------------------------------------
template<typename Container>
inline
Container
concatContainers(const Container& lhs,
                 const Container& rhs) {
  Container result(lhs);
  return concatContainersInPlace(result, rhs);
}

//------------------------------------------------------------------------------
// Repeat
//------------------------------------------------------------------------------
template<typename Container>
inline
Container
repeatContainer(Container& self,
                const unsigned count) {
  Container result(self);
  return repeatContainerInPlace(result, count);
}

//------------------------------------------------------------------------------
// Test if a container contains the given value.
//------------------------------------------------------------------------------
template<typename Container>
inline
int
containsValue(Container& container,
              const typename Container::value_type& value) {
   return ((std::find(container.begin(), container.end(), value) == container.end()) ?
           0 :
           1);
}

template<typename Vptr>
struct PtrComparator {
   bool operator()(const Vptr lhs, const Vptr rhs) {
      return *lhs == *rhs;
   }
};

template<typename Container>
inline
int
containsPtr(Container& container,
            const typename Container::value_type& valuePtr) {
   return ((std::find(container.begin(), container.end(), valuePtr) == container.end()) ?
           0 :
           1);
//    return ((std::find(container.begin(), container.end(), valuePtr, PtrComparator<typename Container::value_type>()) == container.end()) ?
//            0 :
//            1);
}

//------------------------------------------------------------------------------
// Helpers to extract values from a pair.
//------------------------------------------------------------------------------
template<typename Pair>
inline
typename Pair::first_type
extractFirstPairValue(Pair& self) {
  return self.first;
}

template<typename Pair>
inline
typename Pair::second_type
extractSecondPairValue(Pair& self) {
  return self.second;
}

// //------------------------------------------------------------------------------
// // Iterator class to reflect python's iterator semantics.
// //------------------------------------------------------------------------------
// template<typename Container>
// class
// ContainerIterator {
// public:
//   typedef typename Container::iterator iterator;

//   ContainerIterator():
//     mContainerPtr(0),
//     mValue() {}

//   ContainerIterator(Container& container): 
//     mContainerPtr(&container),
//     mValue(container.begin()) {}

//   ContainerIterator(Container& container,
//                     iterator value): 
//     mContainerPtr(&container),
//     mValue(value) {}

//   ContainerIterator* next() const { 
//     iterator nextItr = mValue + 1;
//     if (mValue < mContainerPtr->end()) {
//       return new ContainerIterator(*mContainerPtr, nextItr);
//     } else {
//       PyErr_SetNone(PyExc_StopIteration);
//       return NULL;
//     }
//   }

//   ContainerIterator* __iter__() const { return this; }

//   Container* mContainerPtr;
//   iterator mValue;
// };

// //------------------------------------------------------------------------------
// // Iterator class to reflect python's iterator semantics.
// //------------------------------------------------------------------------------
// template<typename Container>
// inline
// ContainerIterator<Container>*
// nextIterator(ContainerIterator<Container>* itr) {
//   return itr.next();
// }

}

#endif
