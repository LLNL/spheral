//---------------------------------Spheral++----------------------------------//
// Tensor inner double products of Spheral tensor types.
//
// Created by JMO, Thu Oct 15 16:51:37 PDT 2015
//----------------------------------------------------------------------------//
#ifndef __Spheral_innerDoubleProduct_hh__
#define __Spheral_innerDoubleProduct_hh__

namespace Spheral {

//------------------------------------------------------------------------------
// A^{ij} B^{jk}, inner double product of a tensor (rank 2) with a tensor (rank 2).
// This case is actually covered by the A.doubledot(B) method of the rank 2 tensor 
// class.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
innerDoubleProduct(const typename Dimension::Tensor& A,
                   const typename Dimension::Tensor& B) {
  typename Dimension::Scalar C = 0.0;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      C += A(i,j)*B(j,i);
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Scalar
innerDoubleProduct(const typename Dimension::Tensor& A,
                   const typename Dimension::SymTensor& B) {
  typename Dimension::Scalar C = 0.0;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      C += A(i,j)*B(j,i);
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Scalar
innerDoubleProduct(const typename Dimension::SymTensor& A,
                   const typename Dimension::Tensor& B) {
  typename Dimension::Scalar C = 0.0;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      C += A(i,j)*B(j,i);
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Scalar
innerDoubleProduct(const typename Dimension::SymTensor& A,
                   const typename Dimension::SymTensor& B) {
  typename Dimension::Scalar C = 0.0;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      C += A(i,j)*B(j,i);
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ij} B^{jik}, inner double product of a tensor (rank 2) with a tensor (rank 3).
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::Vector
_innerDoubleProduct(const SecondRankTensor& A,
                    const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::Vector C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        C(k) += A(i,j)*B(j,i,k);
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Vector
innerDoubleProduct(const typename Dimension::Tensor& A,
                   const typename Dimension::ThirdRankTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::Vector
innerDoubleProduct(const typename Dimension::SymTensor& A,
             const typename Dimension::ThirdRankTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijk} B^{kj}, inner double product of a tensor (rank 3) with a tensor (rank 2).
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::Vector
_innerDoubleProduct(const typename Dimension::ThirdRankTensor& A,
                    const SecondRankTensor& B) {
  typename Dimension::Vector C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        C(i) += A(i,j,k)*B(k,j);
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Vector
innerDoubleProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::Tensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::Vector
innerDoubleProduct(const typename Dimension::ThirdRankTensor& A,
             const typename Dimension::SymTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijk} B^{kjm}, Inner double product of two rank 3 tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Tensor
innerDoubleProduct(const typename Dimension::ThirdRankTensor& A,
                   const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::Tensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,m) += A(i,j,k)*B(k,j,m);
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ij} B^{jikm}, Inner double product of a rank 2 tensor with a rank 4 tensor.
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::Tensor
_innerDoubleProduct(const SecondRankTensor& A,
                    const typename Dimension::FourthRankTensor& B) {
  typename Dimension::Tensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(k,m) += A(i,j)*B(j,i,k,m);
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerDoubleProduct(const typename Dimension::Tensor& A,
                   const typename Dimension::FourthRankTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerDoubleProduct(const typename Dimension::SymTensor& A,
                   const typename Dimension::FourthRankTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijkm} B^{mk}, Inner double product of a rank 4 tensor with a rank 2 tensor.
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::Tensor
_innerDoubleProduct(const typename Dimension::FourthRankTensor& A,
                    const SecondRankTensor& B) {
  typename Dimension::Tensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          C(i,j) += A(i,j,k,m)*B(m,k);
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerDoubleProduct(const typename Dimension::FourthRankTensor& A,
                   const typename Dimension::Tensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::Tensor
innerDoubleProduct(const typename Dimension::FourthRankTensor& A,
                   const typename Dimension::SymTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijk} B^{kjmn}, Inner double product of a rank 3 tensor with a rank 4 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerDoubleProduct(const typename Dimension::ThirdRankTensor& A,
                   const typename Dimension::FourthRankTensor& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,m,n) += A(i,j,k)*B(k,j,m,n);
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ijkm} B^{mkn}, Inner double product of a rank 4 tensor with a rank 3 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerDoubleProduct(const typename Dimension::FourthRankTensor& A,
                   const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,n) += A(i,j,k,m)*B(m,k,n);
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ijkm} B^{mknp}, Inner double product of two rank 4 tensors.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerDoubleProduct(const typename Dimension::FourthRankTensor& A,
                   const typename Dimension::FourthRankTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            for (size_t p = 0; p != Dimension::nDim; ++p) {
              C(i,j,n,p) += A(i,j,k,m)*B(m,k,n,p);
            }
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ij} B^{jikmn}, Inner double product of a rank 2 tensor with a rank 5 tensor.
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::ThirdRankTensor
_innerDoubleProduct(const SecondRankTensor& A,
                    const typename Dimension::FifthRankTensor& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(k,m,n) += A(i,j)*B(j,i,k,m,n);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerDoubleProduct(const typename Dimension::Tensor& A,
                   const typename Dimension::FifthRankTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerDoubleProduct(const typename Dimension::SymTensor& A,
                   const typename Dimension::FifthRankTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijkmn} B^{nm}, Inner double product of a rank 5 tensor with a rank 2 tensor.
//------------------------------------------------------------------------------
template<typename Dimension, typename SecondRankTensor>
inline
typename Dimension::ThirdRankTensor
_innerDoubleProduct(const typename Dimension::FifthRankTensor& A,
                    const SecondRankTensor& B) {
  typename Dimension::ThirdRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            C(i,j,k) += A(i,j,k,m,n)*B(n,m);
          }
        }
      }
    }
  }
  return C;
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerDoubleProduct(const typename Dimension::FifthRankTensor& A,
                   const typename Dimension::Tensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::Tensor>(A, B);
}

template<typename Dimension>
inline
typename Dimension::ThirdRankTensor
innerDoubleProduct(const typename Dimension::FifthRankTensor& A,
                   const typename Dimension::SymTensor& B) {
  return _innerDoubleProduct<Dimension, typename Dimension::SymTensor>(A, B);
}

//------------------------------------------------------------------------------
// A^{ijk} B^{kjmnp}, Inner double product of a rank 3 tensor with a rank 5 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerDoubleProduct(const typename Dimension::ThirdRankTensor& A,
                   const typename Dimension::FifthRankTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            for (size_t p = 0; p != Dimension::nDim; ++p) {
              C(i,m,n,p) += A(i,j,k)*B(k,j,m,n,p);
            }
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ijkmn} B^{nmp}, Inner double product of a rank 5 tensor with a rank 3 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FourthRankTensor
innerDoubleProduct(const typename Dimension::FifthRankTensor& A,
                   const typename Dimension::ThirdRankTensor& B) {
  typename Dimension::FourthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            for (size_t p = 0; p != Dimension::nDim; ++p) {
              C(i,j,k,p) += A(i,j,k,m,n)*B(n,m,p);
            }
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ijkm} B^{mknpq}, Inner double product of a rank 4 tensor with a rank 5 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FifthRankTensor
innerDoubleProduct(const typename Dimension::FourthRankTensor& A,
                   const typename Dimension::FifthRankTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            for (size_t p = 0; p != Dimension::nDim; ++p) {
              for (size_t q = 0; q != Dimension::nDim; ++q) {
                C(i,j,n,p,q) += A(i,j,k,m)*B(m,k,n,p,q);
              }
            }
          }
        }
      }
    }
  }
  return C;
}

//------------------------------------------------------------------------------
// A^{ijkmn} B^{nmpq}, Inner double product of a rank 5 tensor with a rank 4 tensor.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::FifthRankTensor
innerDoubleProduct(const typename Dimension::FifthRankTensor& A,
                   const typename Dimension::FourthRankTensor& B) {
  typename Dimension::FifthRankTensor C;
  for (size_t i = 0; i != Dimension::nDim; ++i) {
    for (size_t j = 0; j != Dimension::nDim; ++j) {
      for (size_t k = 0; k != Dimension::nDim; ++k) {
        for (size_t m = 0; m != Dimension::nDim; ++m) {
          for (size_t n = 0; n != Dimension::nDim; ++n) {
            for (size_t p = 0; p != Dimension::nDim; ++p) {
              for (size_t q = 0; q != Dimension::nDim; ++q) {
                C(i,j,k,p,q) += A(i,j,k,m,n)*B(n,m,p,q);
              }
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
