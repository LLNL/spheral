//---------------------------------Spheral++----------------------------------//
// Tensor outer products of Spheral vector and tensor types.
//
// Created by JMO, Fri Jul 18 16:42:12 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_outerProduct_hh__
#define __Spheral_outerProduct_hh__

namespace Spheral {

//------------------------------------------------------------------------------
// Outer product with a scalar, just return the scalar times the elements.
//------------------------------------------------------------------------------
template<typename Dimension> inline typename Dimension::Scalar           outerProduct(const typename Dimension::Scalar& A, const typename Dimension::Scalar& B)           { return A*B; }
template<typename Dimension> inline typename Dimension::Vector           outerProduct(const typename Dimension::Scalar& A, const typename Dimension::Vector& B)           { return A*B; }
template<typename Dimension> inline typename Dimension::Tensor           outerProduct(const typename Dimension::Scalar& A, const typename Dimension::Tensor& B)           { return A*B; }
template<typename Dimension> inline typename Dimension::SymTensor        outerProduct(const typename Dimension::Scalar& A, const typename Dimension::SymTensor& B)        { return A*B; }
template<typename Dimension> inline typename Dimension::ThirdRankTensor  outerProduct(const typename Dimension::Scalar& A, const typename Dimension::ThirdRankTensor& B)  { return A*B; }
template<typename Dimension> inline typename Dimension::FourthRankTensor outerProduct(const typename Dimension::Scalar& A, const typename Dimension::FourthRankTensor& B) { return A*B; }
template<typename Dimension> inline typename Dimension::FifthRankTensor  outerProduct(const typename Dimension::Scalar& A, const typename Dimension::FifthRankTensor& B)  { return A*B; }

template<typename Dimension> inline typename Dimension::Vector           outerProduct(const typename Dimension::Vector& A,           const typename Dimension::Scalar& B) { return A*B; }
template<typename Dimension> inline typename Dimension::Tensor           outerProduct(const typename Dimension::Tensor& A,           const typename Dimension::Scalar& B) { return A*B; }
template<typename Dimension> inline typename Dimension::SymTensor        outerProduct(const typename Dimension::SymTensor& A,        const typename Dimension::Scalar& B) { return A*B; }
template<typename Dimension> inline typename Dimension::ThirdRankTensor  outerProduct(const typename Dimension::ThirdRankTensor& A,  const typename Dimension::Scalar& B) { return A*B; }
template<typename Dimension> inline typename Dimension::FourthRankTensor outerProduct(const typename Dimension::FourthRankTensor& A, const typename Dimension::Scalar& B) { return A*B; }
template<typename Dimension> inline typename Dimension::FifthRankTensor  outerProduct(const typename Dimension::FifthRankTensor& A,  const typename Dimension::Scalar& B) { return A*B; }

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

//------------------------------------------------------------------------------
// Outer product of a rank 3 tensor with a Vector.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FourthRankTensor
outerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::Vector& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k,m) = A(i,j,k)*B(m);
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FourthRankTensor
outerProduct(const typename Dimension::Vector& A,
             const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k,m) = A(i)*B(j,k,m);
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// Outer product of two rank 2 tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FourthRankTensor
outerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::Tensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k,m) = A(i,j)*B(k,m);
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FourthRankTensor
outerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::SymTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k,m) = A(i,j)*B(k,m);
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FourthRankTensor
outerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::SymTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k,m) = A(i,j)*B(k,m);
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FourthRankTensor
outerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::Tensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k,m) = A(i,j)*B(k,m);
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// Outer product of a rank 4 tensor with a Vector.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FifthRankTensor
outerProduct(const typename Dimension::FourthRankTensor& A,
             const typename Dimension::Vector& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,m,n) = A(i,j,k,m)*B(n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FifthRankTensor
outerProduct(const typename Dimension::Vector& A,
             const typename Dimension::FourthRankTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,m,n) = A(i)*B(j,k,m,n);
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// Outer product of a rank 3 and rank 2 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FifthRankTensor
outerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::Tensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,m,n) = A(i,j,k)*B(m,n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FifthRankTensor
outerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::SymTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,m,n) = A(i,j,k)*B(m,n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FifthRankTensor
outerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,m,n) = A(i,j)*B(k,m,n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FifthRankTensor
outerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,m,n) = A(i,j)*B(k,m,n);
          }
        }
      }
    }
  }
  return C;
}

}

#endif
