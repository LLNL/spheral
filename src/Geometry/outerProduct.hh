//---------------------------------Spheral++----------------------------------//
// Tensor outer products of Spheral vector and tensor types.
//
// Created by JMO, Fri Jul 18 16:42:12 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_outerProduct_hh__
#define __Spheral_outerProduct_hh__

namespace Geometry {

//------------------------------------------------------------------------------
// Outer product with a scalar, just return the scalar times the elements.
//------------------------------------------------------------------------------
template<typename ValueType>
inline
ValueType
outerProduct(const double& A, const ValueType& B) {
  return A*B;
}

template<typename ValueType>
inline
ValueType
outerProduct(const ValueType& A, const double& B) {
  return A*B;
}

//------------------------------------------------------------------------------
// Vector (rank 1 tensor) outer product, returns a rank 2 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Tensor
outerProduct(const typename Dimension::Vector& A,
             const typename Dimension::Vector& B) {
  return A.dyad(B);
}

//------------------------------------------------------------------------------
// A^{ij} B^{i}, outer product of a rank 2 tensor with a vector.  Returns a 
// rank 3 tensor.
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::ThirdRankTensor
_outerProduct(const SecondRankTensor& A,
             const typename Dimension::Vector& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        C(i,j,k) = A(i,j)*B(k);
      }
    }
  }
  return C;
}

template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::ThirdRankTensor
_outerProduct(const typename Dimension::Vector& A,
             const SecondRankTensor& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        C(i,j,k) = A(i)*B(j,k);
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
outerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::Vector& B) {
  return _outerProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
outerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::Vector& B) {
  return _outerProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
outerProduct(const typename Dimension::Vector& A,
             const typename Dimension::Tensor& B) {
  return _outerProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
outerProduct(const typename Dimension::Vector& A,
             const typename Dimension::SymTensor& B) {
  return _outerProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

}

#endif
