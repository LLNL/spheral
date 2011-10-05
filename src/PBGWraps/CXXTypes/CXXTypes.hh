#ifndef __PBGWRAPS_CXXTYPES__
#define __PBGWRAPS_CXXTYPES__

#include <stdexcept>
#include <iostream>
#include <vector>
#include <string>
#include "Geometry/Dimension.hh"

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

//------------------------------------------------------------------------------
// Index into a vector.
//------------------------------------------------------------------------------
template<typename Value>
inline
Value
indexVector(std::vector<Value>& container,
            const size_t index) {
  try {
    return container.at(index);
  } catch (std::out_of_range) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return Value();
  }
}

template<typename Value>
inline
Value*
indexVectorAsPointer(std::vector<Value>& container,
                     const size_t index) {
  try {
    return &(container.at(index));
  } catch (std::out_of_range) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return NULL;
  }
}

template<typename Value>
inline
Value&
indexVectorAsReference(std::vector<Value>& container,
                       const size_t index) {
  try {
    return container.at(index);
  } catch (std::out_of_range) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return container.front();
  }
}

template<typename ValuePtr>
inline
ValuePtr
indexVectorOfPointers(std::vector<ValuePtr>& container,
                      const size_t index) {
  try {
    return container.at(index);
  } catch (std::out_of_range) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
    return NULL;
  }
}

//------------------------------------------------------------------------------
// Assign to a postion in a container.
//------------------------------------------------------------------------------
template<typename Container>
inline
void
assignToPosition(Container& container, 
                 const size_t index,
                 const typename Container::value_type& value) {
  if (index >= container.size()) {
    PyErr_SetString(PyExc_IndexError, "Container index out of range");
  } else {
    container[index] = value;
  }
}

//------------------------------------------------------------------------------
// Append a value to vector of pointers.
//------------------------------------------------------------------------------
template<typename Container>
inline
void
appendToVectorOfPointers(Container& container, 
                         typename Container::value_type value) {
  container.push_back(value);
}

//------------------------------------------------------------------------------
// Support concatenation of vectors as addition.
//------------------------------------------------------------------------------
template<typename Value>
inline
std::vector<Value>
operator+(const std::vector<Value>& lhs,
          const std::vector<Value>& rhs) {
  std::vector<Value> result;
  std::copy(lhs.begin(), lhs.end(), std::back_inserter(result));
  std::copy(rhs.begin(), rhs.end(), std::back_inserter(result));
  VERIFY(result.size() == lhs.size() + rhs.size());
  return result;
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

#endif
