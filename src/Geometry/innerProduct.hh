//---------------------------------Spheral++----------------------------------//
// Tensor inner products of Spheral vector and tensor types.
//
// Created by JMO, Fri Jul 18 16:42:12 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_innerProduct_hh__
#define __Spheral_innerProduct_hh__

namespace Spheral {

//------------------------------------------------------------------------------
// Inner product with a scalar, just return the scalar times the elements.
//------------------------------------------------------------------------------
template<typename ValueType>
inline
ValueType
innerProduct(const double& A, const ValueType& B) {
  return A*B;
}

template<typename ValueType>
inline
ValueType
innerProduct(const ValueType& B, const double& A) {
  return A*B;
}

//------------------------------------------------------------------------------
// Vector (rank 1 tensor) inner product, returns a rank 0 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
double
innerProduct(const typename Dimension::Vector& A,
             const typename Dimension::Vector& B) {
   return A.dot(B);
}

//------------------------------------------------------------------------------
// A^{ij} B^{j}, inner product of a rank 2 tensor with a vector.  Returns a 
// vector (rank 1 tensor).
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::Vector
_innerProduct(const SecondRankTensor& A,
              const typename Dimension::Vector& B) {
  typename Dimension::Vector C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      C(i) += A(i,j)*B(j);
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Vector
innerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::Vector& B) {
  return _innerProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::Vector
innerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::Vector& B) {
  return _innerProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{j} B^{j,i}, inner product of vector with a rank 2 tensor.  Returns a 
// vector (rank 1 tensor).
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::Vector
_innerProduct(const typename Dimension::Vector& A,
              const SecondRankTensor& B) {
   typename Dimension::Vector C;
   for (size_t i = 0; i != Dimension::nDim; ++i) {
     for (size_t j = 0; j != Dimension::nDim; ++j) {
       C(i) += A(j)*B(j,i);
     }
   }
   return C;
}

template<typename Dimension>
inline
typename Dimension::Vector
innerProduct(const typename Dimension::Vector& A,
             const typename Dimension::Tensor& B) {
  return _innerProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::Vector
innerProduct(const typename Dimension::Vector& A,
             const typename Dimension::SymTensor& B) {
  return _innerProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijk} B^{k}, inner product of a third rank tensor with a vector.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Tensor
innerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::Vector& B) {
  typename Dimension::Tensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        C(i,j) += A(i,j,k)*B(k);
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerProduct(const typename Dimension::Vector& A,
             const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::Tensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        C(i,j) += A(k)*B(k,i,j);
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ij} B^{jk}, inner product of a tensor (rank 2) with a tensor (rank 2).
// This case is actually covered by the A.dot(B) method of the rank 2 tensor 
// class.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Tensor
innerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::Tensor& B) {
   return A.dot(B);
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::SymTensor& B) {
   return A.dot(B);
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::SymTensor& B) {
   return A.dot(B);
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::Tensor& B) {
   return A.dot(B);
}

//------------------------------------------------------------------------------
// A^{ij} B^{jkl}, inner product of a tensor (rank 2) with a tensor (rank 3).
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::ThirdRankTensor
_innerProduct(const SecondRankTensor& A,
              const typename Dimension::ThirdRankTensor& B) {
   typename Dimension::ThirdRankTensor C;
   for (size_t i = 0; i != Dimension::nDim; ++i) {
     for (size_t k = 0; k != Dimension::nDim; ++k) {
       for (size_t l = 0; l != Dimension::nDim; ++l) {
         for (size_t j = 0; j != Dimension::nDim; ++j) {
           C(i,k,l) += A(i,j)*B(j,k,l);
         }
       }
     }
   }
   return C;
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::ThirdRankTensor& B) {
  return _innerProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::ThirdRankTensor& B) {
  return _innerProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijk} B^{kl}, inner product of a tensor (rank 3) with a tensor (rank 2).
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::ThirdRankTensor
_innerProduct(const typename Dimension::ThirdRankTensor& A,
              const SecondRankTensor& B) {
   typename Dimension::ThirdRankTensor C;
   for (size_t i = 0; i != Dimension::nDim; ++i) {
     for (size_t j = 0; j != Dimension::nDim; ++j) {
       for (size_t l = 0; l != Dimension::nDim; ++l) {
         for (size_t k = 0; k != Dimension::nDim; ++k) {
           C(i,j,l) += A(i,j,k)*B(k,l);
         }
       }
     }
   }
   return C;
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::Tensor& B) {
  return _innerProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::SymTensor& B) {
  return _innerProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// Inner product of a rank 3 tensor with a rank 3 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,m,n) += A(i,j,k)*B(k,m,n);
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// Inner product of a rank 4 tensor with a Vector.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerProduct(const typename Dimension::FourthRankTensor& A,
             const typename Dimension::Vector& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k) += A(i,j,k,m)*B(m);
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerProduct(const typename Dimension::Vector& A,
             const typename Dimension::FourthRankTensor& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j,k) += A(m)*B(m,i,j,k);
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// Inner product of a rank 4 tensor with a rank 2 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerProduct(const typename Dimension::FourthRankTensor& A,
             const typename Dimension::Tensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,n) += A(i,j,k,m)*B(m,n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerProduct(const typename Dimension::FourthRankTensor& A,
             const typename Dimension::SymTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k,n) += A(i,j,k,m)*B(m,n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerProduct(const typename Dimension::Tensor& A,
             const typename Dimension::FourthRankTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,k,m,n) += A(i,j)*B(j,k,m,n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::FourthRankTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,k,m,n) += A(i,j)*B(j,k,m,n);
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// Inner product of a rank 4 tensor with a rank 3 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FifthRankTensor
innerProduct(const typename Dimension::FourthRankTensor& A,
             const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            for (size_t p = 0; p != Dimension::nDim; ++p) {
              C(i,j,k,n,p) += A(i,j,k,m)*B(m,n,p);
            }
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
innerProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::FourthRankTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            for (size_t p = 0; p != Dimension::nDim; ++p) {
              C(i,j,m,n,p) += A(i,j,k)*B(k,m,n,p);
            }
          }
        }
      }
    }
  }
  return C;
}

}
#endif
