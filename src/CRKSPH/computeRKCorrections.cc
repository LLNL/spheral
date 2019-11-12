//------------------------------------------------------------------------------
// Compute the RK correction terms
//------------------------------------------------------------------------------

#include <array>
#include "computeRKCorrections.hh"
#include "Eigen/Dense"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {

// template<int nDim>
// class GeomSixthRankTensor : public RankNTensor<nDim, 6, SixthRankTensor<nDim> {
// public:
//   typedef typename RankNTensor<nDim, 6, GeomSixthRankTensor>::size_type size_type;
//   static const size_type numElements;
  
//   // Useful static member data.
//   static const GeomFourthRankTensor zero;
  
//   // Constructor
//   SixthRankTensor():
//     RankNTensor<nDim, 5, GeomFifthRankTensor>() {
//   }
//   SixthRankTensor(const double val):
//     RankNTensor<nDim, 5, GeomFifthRankTensor>(val) {
//   }
  
//   double operator()(const size_type i, const size_type j, const size_type l, const size_type m, const size_type n, const size_type o) {
//     REQUIRE(i < nDim and j < nDim and k < nDim and m < nDim and n < nDim and o < nDim);
//     return mElements[((((i*nDim + j)*nDim + k)*nDim + m)*nDim + n)*nDim + o];
//   }
  
//   double& operator()(const size_type i, const size_type j, const size_type l, const size_type m, const size_type n, const size_type o) {
//     REQUIRE(i < nDim and j < nDim and k < nDim and m < nDim and n < nDim and o < nDim);
//     return mElements[((((i*nDim + j)*nDim + k)*nDim + m)*nDim + n)*nDim + o];
//   }
// }
  
// template<> const unsigned GeomSixthRankTensor<1>::numElements = 1;
// template<> const GeomSixthRankTensor<1> GeomSixthRankTensor<1>::zero = GeomSixthRankTensor<1>(0.0);
// template<> const unsigned GeomSixthRankTensor<2>::numElements = 64;
// template<> const GeomSixthRankTensor<2> GeomSixthRankTensor<2>::zero = GeomSixthRankTensor<2>(0.0);
// template<> const unsigned GeomSixthRankTensor<3>::numElements = 729;
// template<> const GeomSixthRankTensor<3> GeomSixthRankTensor<3>::zero = GeomSixthRankTensor<3>(0.0);
    

// template<typename Dimension>
// using SixthRankTensor = 
// template<typename Dimension>
// using SeventhRankTensor = std::array<double, Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim>;
// template<typename Dimension>
// using EighthRankTensor = std::array<double, Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim * Dimension::nDim>;

//------------------------------------------------------------------------------
// Struct to carry around the moment values
//------------------------------------------------------------------------------
template<typename Dimension>
struct RKMomentValues {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor Tensor3;
  typedef typename Dimension::FourthRankTensor Tensor4;
  typedef typename Dimension::FifthRankTensor Tensor5;
  typedef SixthRankTensor<Dimension> Tensor6;
  typedef SeventhRankTensor<Dimension> Tensor7;
  typedef EighthRankTensor<Dimension> Tensor8;

  // Moments
  Scalar m0;
  Vector m1;
  Tensor m2;
  Tensor3 m3;
  Tensor4 m4;
  Tensor5 m5;
  Tensor6 m6;

  // Gradients
  Vector dm0;
  Tensor dm1;
  Tensor3 dm2;
  Tensor4 dm3;
  Tensor5 dm4;
  Tensor6 dm5;
  Tensor7 dm6;

  // Hessians
  Tensor ddm0;
  Tensor3 ddm1;
  Tensor4 ddm2;
  Tensor5 ddm3;
  Tensor6 ddm4;
  Tensor7 ddm5;
  Tensor8 ddm6;
};

//------------------------------------------------------------------------------
// Zeroth order moments
//------------------------------------------------------------------------------
template<typename Dimension>
void
addToMomentsZeroth(const typename Dimension::Vector& eta,
                   const typename Dimension::Tensor& deta,
                   const typename Dimension::ThirdRankTensor& ddeta,
                   const typename Dimension::Scalar& w,
                   const typename Dimension::Vector& dw,
                   const typename Dimension::Tensor& ddw,
                   const typename Dimension::Scalar& v,
                   const CRKOrder correctionOrder,
                   RKMomentValues<Dimension>& moments) {
  constexpr auto dim = Dimension::nDim;

  // Moments
  moments.m0 += v*w;
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    moments.dm0[k1] += v*dw[k1];
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      moments.ddm0[k1][k2] += v*ddw[k1][k2];
    }
  }
}

//------------------------------------------------------------------------------
// Linear moments
//------------------------------------------------------------------------------
template<typename Dimension>
void
addToMomentsLinear(const typename Dimension::Vector& eta,
                   const typename Dimension::Tensor& deta,
                   const typename Dimension::ThirdRankTensor& ddeta,
                   const typename Dimension::Scalar& w,
                   const typename Dimension::Vector& dw,
                   const typename Dimension::Tensor& ddw,
                   const typename Dimension::Scalar& v,
                   const CRKOrder correctionOrder,
                   RKMomentValues<Dimension>& moments) {
  constexpr auto dim = Dimension::nDim;

  // Moments
  moments.m0 += v*w;
  for (auto q1 = 0; q1 < dim; ++q1) {
    moments.m1[q1] += v*w*eta[q1];
    for (auto q2 = 0; q2 < dim; ++q2) {
      moments.m2[q1][q2] += v*w*eta[q1]*eta[q2];
    }
  }
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    moments.dm0[k1] += v*dw[k1];
    for (auto q1 = 0; q1 < dim; ++q1) {
      moments.dm1[q1][k1] += v*dw[k1]*eta[q1] + v*w*deta[q1][k1];
      for (auto q2 = 0; q2 < dim; ++q2) {
        moments.dm2[q1][q2][k1] += v*dw[k1]*eta[q1]*eta[q2] + v*w*eta[q2]*deta[q1][k1] + v*w*eta[q1]*deta[q2][k1];
      }
    }
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      moments.ddm0[k1][k2] += v*ddw[k1][k2];
      for (auto q1 = 0; q1 < dim; ++q1) {
        moments.ddm1[q1][k1][k2] += v*ddw[k1][k2]*eta[q1] + v*w*ddeta[q1][k1][k2] + v*dw[k2]*deta[q1][k1] + v*dw[k1]*deta[q1][k2];
        for (auto q2 = 0; q2 < dim; ++q2) {
          moments.ddm2[q1][q2][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2] + v*w*eta[q2]*ddeta[q1][k1][k2] + v*w*eta[q1]*ddeta[q2][k1][k2] + v*dw[k2]*eta[q2]*deta[q1][k1] + v*dw[k1]*eta[q2]*deta[q1][k2] + v*dw[k2]*eta[q1]*deta[q2][k1] + v*w*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*deta[q2][k2] + v*w*deta[q1][k1]*deta[q2][k2];
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Quadratic moments
//------------------------------------------------------------------------------
template<typename Dimension>
void
addToMomentsQuadratic(const typename Dimension::Vector& eta,
                   const typename Dimension::Tensor& deta,
                   const typename Dimension::ThirdRankTensor& ddeta,
                   const typename Dimension::Scalar& w,
                   const typename Dimension::Vector& dw,
                   const typename Dimension::Tensor& ddw,
                   const typename Dimension::Scalar& v,
                   const CRKOrder correctionOrder,
                   RKMomentValues<Dimension>& moments) {
  constexpr auto dim = Dimension::nDim;
  
  // Moments
  moments.m0 += v*w;
  for (auto q1 = 0; q1 < dim; ++q1) {
    moments.m1[q1] += v*w*eta[q1];
    for (auto q2 = 0; q2 < dim; ++q2) {
      moments.m2[q1][q2] += v*w*eta[q1]*eta[q2];
      for (auto q3 = 0; q3 < dim; ++q3) {
        moments.m3[q1][q2][q3] += v*w*eta[q1]*eta[q2]*eta[q3];
        for (auto q4 = 0; q4 < dim; ++q4) {
          moments.m4[q1][q2][q3][q4] += v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4];
        }
      }
    }
  }
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    moments.dm0[k1] += v*dw[k1];
    for (auto q1 = 0; q1 < dim; ++q1) {
      moments.dm1[q1][k1] += v*dw[k1]*eta[q1] + v*w*deta[q1][k1];
      for (auto q2 = 0; q2 < dim; ++q2) {
        moments.dm2[q1][q2][k1] += v*dw[k1]*eta[q1]*eta[q2] + v*w*eta[q2]*deta[q1][k1] + v*w*eta[q1]*deta[q2][k1];
        for (auto q3 = 0; q3 < dim; ++q3) {
          moments.dm3[q1][q2][q3][k1] += v*dw[k1]*eta[q1]*eta[q2]*eta[q3] + v*w*eta[q2]*eta[q3]*deta[q1][k1] + v*w*eta[q1]*eta[q3]*deta[q2][k1] + v*w*eta[q1]*eta[q2]*deta[q3][k1];
          for (auto q4 = 0; q4 < dim; ++q4) {
            moments.dm4[q1][q2][q3][q4][k1] += v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q4] + v*w*eta[q2]*eta[q3]*eta[q4]*deta[q1][k1] + v*w*eta[q1]*eta[q3]*eta[q4]*deta[q2][k1] + v*w*eta[q1]*eta[q2]*eta[q4]*deta[q3][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*deta[q4][k1];
          }
        }
      }
    }
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      moments.ddm0[k1][k2] += v*ddw[k1][k2];
      for (auto q1 = 0; q1 < dim; ++q1) {
        moments.ddm1[q1][k1][k2] += v*ddw[k1][k2]*eta[q1] + v*w*ddeta[q1][k1][k2] + v*dw[k2]*deta[q1][k1] + v*dw[k1]*deta[q1][k2];
        for (auto q2 = 0; q2 < dim; ++q2) {
          moments.ddm2[q1][q2][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2] + v*w*eta[q2]*ddeta[q1][k1][k2] + v*w*eta[q1]*ddeta[q2][k1][k2] + v*dw[k2]*eta[q2]*deta[q1][k1] + v*dw[k1]*eta[q2]*deta[q1][k2] + v*dw[k2]*eta[q1]*deta[q2][k1] + v*w*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*deta[q2][k2] + v*w*deta[q1][k1]*deta[q2][k2];
          for (auto q3 = 0; q3 < dim; ++q3) {
            moments.ddm3[q1][q2][q3][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2]*eta[q3] + v*w*eta[q2]*eta[q3]*ddeta[q1][k1][k2] + v*w*eta[q1]*eta[q3]*ddeta[q2][k1][k2] + v*w*eta[q1]*eta[q2]*ddeta[q3][k1][k2] + v*dw[k2]*eta[q2]*eta[q3]*deta[q1][k1] + v*dw[k1]*eta[q2]*eta[q3]*deta[q1][k2] + v*dw[k2]*eta[q1]*eta[q3]*deta[q2][k1] + v*w*eta[q3]*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*eta[q3]*deta[q2][k2] + v*w*eta[q3]*deta[q1][k1]*deta[q2][k2] + v*dw[k2]*eta[q1]*eta[q2]*deta[q3][k1] + v*w*eta[q2]*deta[q1][k2]*deta[q3][k1] + v*w*eta[q1]*deta[q2][k2]*deta[q3][k1] + v*dw[k1]*eta[q1]*eta[q2]*deta[q3][k2] + v*w*eta[q2]*deta[q1][k1]*deta[q3][k2] + v*w*eta[q1]*deta[q2][k1]*deta[q3][k2];
            for (auto q4 = 0; q4 < dim; ++q4) {
              moments.ddm4[q1][q2][q3][q4][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2]*eta[q3]*eta[q4] + v*w*eta[q2]*eta[q3]*eta[q4]*ddeta[q1][k1][k2] + v*w*eta[q1]*eta[q3]*eta[q4]*ddeta[q2][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q4]*ddeta[q3][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*ddeta[q4][k1][k2] + v*dw[k2]*eta[q2]*eta[q3]*eta[q4]*deta[q1][k1] + v*dw[k1]*eta[q2]*eta[q3]*eta[q4]*deta[q1][k2] + v*dw[k2]*eta[q1]*eta[q3]*eta[q4]*deta[q2][k1] + v*w*eta[q3]*eta[q4]*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*eta[q3]*eta[q4]*deta[q2][k2] + v*w*eta[q3]*eta[q4]*deta[q1][k1]*deta[q2][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q4]*deta[q3][k1] + v*w*eta[q2]*eta[q4]*deta[q1][k2]*deta[q3][k1] + v*w*eta[q1]*eta[q4]*deta[q2][k2]*deta[q3][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q4]*deta[q3][k2] + v*w*eta[q2]*eta[q4]*deta[q1][k1]*deta[q3][k2] + v*w*eta[q1]*eta[q4]*deta[q2][k1]*deta[q3][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q3]*deta[q4][k1] + v*w*eta[q2]*eta[q3]*deta[q1][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q3]*deta[q2][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q2]*deta[q3][k2]*deta[q4][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*deta[q4][k2] + v*w*eta[q2]*eta[q3]*deta[q1][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q3]*deta[q2][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q2]*deta[q3][k1]*deta[q4][k2];
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Cubic moments
//------------------------------------------------------------------------------
template<typename Dimension>
void
addToMomentsCubic(const typename Dimension::Vector& eta,
                  const typename Dimension::Tensor& deta,
                  const typename Dimension::ThirdRankTensor& ddeta,
                  const typename Dimension::Scalar& w,
                  const typename Dimension::Vector& dw,
                  const typename Dimension::Tensor& ddw,
                  const typename Dimension::Scalar& v,
                  const CRKOrder correctionOrder,
                  RKMomentValues<Dimension>& moments) {
  constexpr auto dim = Dimension::nDim;

  // Moments
  moments.m0 += v*w;
  for (auto q1 = 0; q1 < dim; ++q1) {
    moments.m1[q1] += v*w*eta[q1];
    for (auto q2 = 0; q2 < dim; ++q2) {
      moments.m2[q1][q2] += v*w*eta[q1]*eta[q2];
      for (auto q3 = 0; q3 < dim; ++q3) {
        moments.m3[q1][q2][q3] += v*w*eta[q1]*eta[q2]*eta[q3];
        for (auto q4 = 0; q4 < dim; ++q4) {
          moments.m4[q1][q2][q3][q4] += v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4];
          for (auto q5 = 0; q5 < dim; ++q5) {
            moments.m5[q1][q2][q3][q4][q5] += v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5];
            for (auto q6 = 0; q6 < dim; ++q6) {
              moments.m6[q1][q2][q3][q4][q5][q6] += v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*eta[q6];
            }
          }
        }
      }
    }
  }
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    moments.dm0[k1] += v*dw[k1];
    for (auto q1 = 0; q1 < dim; ++q1) {
      moments.dm1[q1][k1] += v*dw[k1]*eta[q1] + v*w*deta[q1][k1];
      for (auto q2 = 0; q2 < dim; ++q2) {
        moments.dm2[q1][q2][k1] += v*dw[k1]*eta[q1]*eta[q2] + v*w*eta[q2]*deta[q1][k1] + v*w*eta[q1]*deta[q2][k1];
        for (auto q3 = 0; q3 < dim; ++q3) {
          moments.dm3[q1][q2][q3][k1] += v*dw[k1]*eta[q1]*eta[q2]*eta[q3] + v*w*eta[q2]*eta[q3]*deta[q1][k1] + v*w*eta[q1]*eta[q3]*deta[q2][k1] + v*w*eta[q1]*eta[q2]*deta[q3][k1];
          for (auto q4 = 0; q4 < dim; ++q4) {
            moments.dm4[q1][q2][q3][q4][k1] += v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q4] + v*w*eta[q2]*eta[q3]*eta[q4]*deta[q1][k1] + v*w*eta[q1]*eta[q3]*eta[q4]*deta[q2][k1] + v*w*eta[q1]*eta[q2]*eta[q4]*deta[q3][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*deta[q4][k1];
            for (auto q5 = 0; q5 < dim; ++q5) {
              moments.dm5[q1][q2][q3][q4][q5][k1] += v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q1][k1] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q5]*deta[q2][k1] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q5]*deta[q3][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q5]*deta[q4][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*deta[q5][k1];
              for (auto q6 = 0; q6 < dim; ++q6) {
                moments.dm6[q1][q2][q3][q4][q5][q6][k1] += v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*eta[q6] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q1][k1] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q2][k1] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q5]*eta[q6]*deta[q3][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q5]*eta[q6]*deta[q4][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q6]*deta[q5][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q6][k1];
              }
            }
          }
        }
      }
    }
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      moments.ddm0[k1][k2] += v*ddw[k1][k2];
      for (auto q1 = 0; q1 < dim; ++q1) {
        moments.ddm1[q1][k1][k2] += v*ddw[k1][k2]*eta[q1] + v*w*ddeta[q1][k1][k2] + v*dw[k2]*deta[q1][k1] + v*dw[k1]*deta[q1][k2];
        for (auto q2 = 0; q2 < dim; ++q2) {
          moments.ddm2[q1][q2][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2] + v*w*eta[q2]*ddeta[q1][k1][k2] + v*w*eta[q1]*ddeta[q2][k1][k2] + v*dw[k2]*eta[q2]*deta[q1][k1] + v*dw[k1]*eta[q2]*deta[q1][k2] + v*dw[k2]*eta[q1]*deta[q2][k1] + v*w*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*deta[q2][k2] + v*w*deta[q1][k1]*deta[q2][k2];
          for (auto q3 = 0; q3 < dim; ++q3) {
            moments.ddm3[q1][q2][q3][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2]*eta[q3] + v*w*eta[q2]*eta[q3]*ddeta[q1][k1][k2] + v*w*eta[q1]*eta[q3]*ddeta[q2][k1][k2] + v*w*eta[q1]*eta[q2]*ddeta[q3][k1][k2] + v*dw[k2]*eta[q2]*eta[q3]*deta[q1][k1] + v*dw[k1]*eta[q2]*eta[q3]*deta[q1][k2] + v*dw[k2]*eta[q1]*eta[q3]*deta[q2][k1] + v*w*eta[q3]*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*eta[q3]*deta[q2][k2] + v*w*eta[q3]*deta[q1][k1]*deta[q2][k2] + v*dw[k2]*eta[q1]*eta[q2]*deta[q3][k1] + v*w*eta[q2]*deta[q1][k2]*deta[q3][k1] + v*w*eta[q1]*deta[q2][k2]*deta[q3][k1] + v*dw[k1]*eta[q1]*eta[q2]*deta[q3][k2] + v*w*eta[q2]*deta[q1][k1]*deta[q3][k2] + v*w*eta[q1]*deta[q2][k1]*deta[q3][k2];
            for (auto q4 = 0; q4 < dim; ++q4) {
              moments.ddm4[q1][q2][q3][q4][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2]*eta[q3]*eta[q4] + v*w*eta[q2]*eta[q3]*eta[q4]*ddeta[q1][k1][k2] + v*w*eta[q1]*eta[q3]*eta[q4]*ddeta[q2][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q4]*ddeta[q3][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*ddeta[q4][k1][k2] + v*dw[k2]*eta[q2]*eta[q3]*eta[q4]*deta[q1][k1] + v*dw[k1]*eta[q2]*eta[q3]*eta[q4]*deta[q1][k2] + v*dw[k2]*eta[q1]*eta[q3]*eta[q4]*deta[q2][k1] + v*w*eta[q3]*eta[q4]*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*eta[q3]*eta[q4]*deta[q2][k2] + v*w*eta[q3]*eta[q4]*deta[q1][k1]*deta[q2][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q4]*deta[q3][k1] + v*w*eta[q2]*eta[q4]*deta[q1][k2]*deta[q3][k1] + v*w*eta[q1]*eta[q4]*deta[q2][k2]*deta[q3][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q4]*deta[q3][k2] + v*w*eta[q2]*eta[q4]*deta[q1][k1]*deta[q3][k2] + v*w*eta[q1]*eta[q4]*deta[q2][k1]*deta[q3][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q3]*deta[q4][k1] + v*w*eta[q2]*eta[q3]*deta[q1][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q3]*deta[q2][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q2]*deta[q3][k2]*deta[q4][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*deta[q4][k2] + v*w*eta[q2]*eta[q3]*deta[q1][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q3]*deta[q2][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q2]*deta[q3][k1]*deta[q4][k2];
              for (auto q5 = 0; q5 < dim; ++q5) {
                moments.ddm5[q1][q2][q3][q4][q5][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q5]*ddeta[q1][k1][k2] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q5]*ddeta[q2][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q5]*ddeta[q3][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q5]*ddeta[q4][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*ddeta[q5][k1][k2] + v*dw[k2]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q1][k1] + v*dw[k1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q1][k2] + v*dw[k2]*eta[q1]*eta[q3]*eta[q4]*eta[q5]*deta[q2][k1] + v*w*eta[q3]*eta[q4]*eta[q5]*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*eta[q3]*eta[q4]*eta[q5]*deta[q2][k2] + v*w*eta[q3]*eta[q4]*eta[q5]*deta[q1][k1]*deta[q2][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q4]*eta[q5]*deta[q3][k1] + v*w*eta[q2]*eta[q4]*eta[q5]*deta[q1][k2]*deta[q3][k1] + v*w*eta[q1]*eta[q4]*eta[q5]*deta[q2][k2]*deta[q3][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q4]*eta[q5]*deta[q3][k2] + v*w*eta[q2]*eta[q4]*eta[q5]*deta[q1][k1]*deta[q3][k2] + v*w*eta[q1]*eta[q4]*eta[q5]*deta[q2][k1]*deta[q3][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q3]*eta[q5]*deta[q4][k1] + v*w*eta[q2]*eta[q3]*eta[q5]*deta[q1][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q3]*eta[q5]*deta[q2][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q2]*eta[q5]*deta[q3][k2]*deta[q4][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q5]*deta[q4][k2] + v*w*eta[q2]*eta[q3]*eta[q5]*deta[q1][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q3]*eta[q5]*deta[q2][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q2]*eta[q5]*deta[q3][k1]*deta[q4][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*deta[q5][k1] + v*w*eta[q2]*eta[q3]*eta[q4]*deta[q1][k2]*deta[q5][k1] + v*w*eta[q1]*eta[q3]*eta[q4]*deta[q2][k2]*deta[q5][k1] + v*w*eta[q1]*eta[q2]*eta[q4]*deta[q3][k2]*deta[q5][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*deta[q4][k2]*deta[q5][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*deta[q5][k2] + v*w*eta[q2]*eta[q3]*eta[q4]*deta[q1][k1]*deta[q5][k2] + v*w*eta[q1]*eta[q3]*eta[q4]*deta[q2][k1]*deta[q5][k2] + v*w*eta[q1]*eta[q2]*eta[q4]*deta[q3][k1]*deta[q5][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*deta[q4][k1]*deta[q5][k2];
                for (auto q6 = 0; q6 < dim; ++q6) {
                  moments.ddm6[q1][q2][q3][q4][q5][q6][k1][k2] += v*ddw[k1][k2]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*eta[q6] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*ddeta[q1][k1][k2] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*ddeta[q2][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q5]*eta[q6]*ddeta[q3][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q5]*eta[q6]*ddeta[q4][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q6]*ddeta[q5][k1][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*ddeta[q6][k1][k2] + v*dw[k2]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q1][k1] + v*dw[k1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q1][k2] + v*dw[k2]*eta[q1]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q2][k1] + v*w*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q1][k2]*deta[q2][k1] + v*dw[k1]*eta[q1]*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q2][k2] + v*w*eta[q3]*eta[q4]*eta[q5]*eta[q6]*deta[q1][k1]*deta[q2][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q4]*eta[q5]*eta[q6]*deta[q3][k1] + v*w*eta[q2]*eta[q4]*eta[q5]*eta[q6]*deta[q1][k2]*deta[q3][k1] + v*w*eta[q1]*eta[q4]*eta[q5]*eta[q6]*deta[q2][k2]*deta[q3][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q4]*eta[q5]*eta[q6]*deta[q3][k2] + v*w*eta[q2]*eta[q4]*eta[q5]*eta[q6]*deta[q1][k1]*deta[q3][k2] + v*w*eta[q1]*eta[q4]*eta[q5]*eta[q6]*deta[q2][k1]*deta[q3][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q3]*eta[q5]*eta[q6]*deta[q4][k1] + v*w*eta[q2]*eta[q3]*eta[q5]*eta[q6]*deta[q1][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q3]*eta[q5]*eta[q6]*deta[q2][k2]*deta[q4][k1] + v*w*eta[q1]*eta[q2]*eta[q5]*eta[q6]*deta[q3][k2]*deta[q4][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q5]*eta[q6]*deta[q4][k2] + v*w*eta[q2]*eta[q3]*eta[q5]*eta[q6]*deta[q1][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q3]*eta[q5]*eta[q6]*deta[q2][k1]*deta[q4][k2] + v*w*eta[q1]*eta[q2]*eta[q5]*eta[q6]*deta[q3][k1]*deta[q4][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q6]*deta[q5][k1] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q6]*deta[q1][k2]*deta[q5][k1] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q6]*deta[q2][k2]*deta[q5][k1] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q6]*deta[q3][k2]*deta[q5][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q6]*deta[q4][k2]*deta[q5][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q6]*deta[q5][k2] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q6]*deta[q1][k1]*deta[q5][k2] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q6]*deta[q2][k1]*deta[q5][k2] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q6]*deta[q3][k1]*deta[q5][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q6]*deta[q4][k1]*deta[q5][k2] + v*dw[k2]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q6][k1] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q1][k2]*deta[q6][k1] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q5]*deta[q2][k2]*deta[q6][k1] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q5]*deta[q3][k2]*deta[q6][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q5]*deta[q4][k2]*deta[q6][k1] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*deta[q5][k2]*deta[q6][k1] + v*dw[k1]*eta[q1]*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q6][k2] + v*w*eta[q2]*eta[q3]*eta[q4]*eta[q5]*deta[q1][k1]*deta[q6][k2] + v*w*eta[q1]*eta[q3]*eta[q4]*eta[q5]*deta[q2][k1]*deta[q6][k2] + v*w*eta[q1]*eta[q2]*eta[q4]*eta[q5]*deta[q3][k1]*deta[q6][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q5]*deta[q4][k1]*deta[q6][k2] + v*w*eta[q1]*eta[q2]*eta[q3]*eta[q4]*deta[q5][k1]*deta[q6][k2];
                }
              }
            }
          }
        }
      }
    }
  }
}

//------------------------------------------------------------------------------
// Zeroth order corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeZerothCorrections(const RKMomentValues<Dimension>& moments,
                         typename Dimension::Scalar& a,
                         typename Dimension::Vector& da,
                         typename Dimension::Tensor& dda) {
  ASSERT2(false, "not implemented");
}

//------------------------------------------------------------------------------
// Linear corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeLinearCorrections(const RKMomentValues<Dimension>& moments,
                         typename Dimension::Scalar& a,
                         typename Dimension::Vector& b,
                         typename Dimension::Vector& da,
                         typename Dimension::Tensor& db,
                         typename Dimension::Tensor& dda,
                         typename Dimension::ThirdRankTensor& ddb) {
  constexpr auto dim = Dimension::nDim;
  
  // Get m matrix
  constexpr auto size = Dimension::nDim + 1;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  
  // Get moment values
  const auto& m0 = moments.m0;
  const auto& m1 = moments.m1;
  const auto& m2 = moments.m2;
  const auto& dm0 = moments.dm0;
  const auto& dm1 = moments.dm1;
  const auto& dm2 = moments.dm2;
  const auto& ddm0 = moments.ddm0;
  const auto& ddm1 = moments.ddm1;
  const auto& ddm2 = moments.ddm2;

  // Initialize Eigen stuff
  MatrixType matrix;
  VectorType rhs;
  
  // Put m values into matrix
  matrix(0, 0) = m0;
  for (auto q1 = 0; q1 < dim; ++q1) {
    matrix(0, q1+1) = m1[q1];
    matrix(q1+1, 0) = m1[q1];
    for (auto q2 = 0; q2 < dim; ++q2) {
      matrix(q1+1,q2+1) = m2[q1][q2];
    }
  }

  // Perform decomposition
  auto solver = matrix.colPivHouseholderQr();

  // Solve for values
  {
    rhs(0) = 1.;
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      rhs(row) = 0.;
    }
    // Solve
    VectorType solution = solver.solve(rhs);

    // Get result
    a = solution(0);
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      b(p1) = solution(row);
    }
  }
  
  // Solve for derivatives
  rhs.setZero();
  for (auto k1 = 0; k1 < dim; ++k1) {
    // Row 0 
    rhs(0) -= a*dm0[k1];
    for (auto q1 = 0; q1 < dim; ++q1) {
      rhs(0) -= b[q1]*dm1[q1][k1];
    }

    // Row 1-4
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      rhs(row) -= a*dm1[p1][k1];
      for (auto q1 = 0; q1 < dim; ++q1) {
        rhs(row) -= b[q1]*dm2[q1][p1][k1];
      }
    }

    // Solve
    VectorType solution = solver.solve(rhs);

    // Get result
    da[k1] = solution(0);
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      db[p1][k1] = solution(row);
    }
  }
  
  // Solve for second derivatives
  rhs.setZero();
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      // Row 0
      rhs(0) -= m0*dda[k1][k2] + a*ddm0[k1][k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2];
      for (auto q1 = 0; q1 < dim; ++q1) {
        rhs(0) -= m1[q1]*ddb[q1][k1][k2] + b[q1]*ddm1[q1][k1][k2] + db[q1][k2]*dm1[q1][k1] + db[q1][k1]*dm1[q1][k2];
      }

      // Row 1-4
      for (auto p1 = 0; p1 < dim; ++p1) {
        const auto row = p1 + 1;
        rhs(row) -= a*ddm1[p1][k1][k2] + da[k2]*dm1[p1][k1] + da[k1]*dm1[p1][k2];
        for (auto q1 = 0; q1 < dim; ++q1) {
          rhs(row) -= b[q1]*ddm2[q1][p1][k1][k2]  + db[q1][k2]*dm2[q1][p1][k1] + db[q1][k1]*dm2[q1][p1][k2];
        }
      }

      // Solve
      VectorType solution = solver.solve(rhs);

      // Get result
      dda[k1][k2] = solution(0);
      for (auto p1 = 0; p1 < dim; ++p1) {
        const auto row = p1 + 1;
        ddb[p1][k1][k2] = solution(p1);
      }
    }
  }
}

//------------------------------------------------------------------------------
// Quadratic corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeQuadraticCorrections(const RKMomentValues<Dimension>& moments,
                            typename Dimension::Scalar& a,
                            typename Dimension::Vector& b,
                            typename Dimension::Tensor& c,
                            typename Dimension::Vector& da,
                            typename Dimension::Tensor& db,
                            typename Dimension::ThirdRankTensor& dc,
                            typename Dimension::Tensor& dda,
                            typename Dimension::ThirdRankTensor& ddb,
                            typename Dimension::FourthRankTensor& ddc) {
  ASSERT2(false, "not implemented");
}

//------------------------------------------------------------------------------
// Cubic corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeCubicCorrections(const RKMomentValues<Dimension>& moments,
                        typename Dimension::Scalar& a,
                        typename Dimension::Vector& b,
                        typename Dimension::Tensor& c,
                        typename Dimension::Tensor& d,
                        typename Dimension::Vector& da,
                        typename Dimension::Tensor& db,
                        typename Dimension::ThirdRankTensor& dc,
                        typename Dimension::FourthRankTensor& dd,
                        typename Dimension::Tensor& dda,
                        typename Dimension::ThirdRankTensor& ddb,
                        typename Dimension::FourthRankTensor& ddc,
                        typename Dimension::FifthRankTensor& ddd) {
  ASSERT2(false, "not implemented");
}

//------------------------------------------------------------------------------
// Compute the corrections
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeRKCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                     const TableKernel<Dimension>& W,
                     const FieldList<Dimension, typename Dimension::Scalar>& volume,
                     const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     const CRKOrder correctionOrder,
                     FieldList<Dimension, typename Dimension::Scalar>& A,
                     FieldList<Dimension, typename Dimension::Vector>& B,
                     FieldList<Dimension, typename Dimension::Tensor>& C,
                     FieldList<Dimension, typename Dimension::Tensor>& D,
                     FieldList<Dimension, typename Dimension::Vector>& gradA,
                     FieldList<Dimension, typename Dimension::Tensor>& gradB,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& gradD,
                     FieldList<Dimension, typename Dimension::Tensor>& hessA,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& hessB,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& hessC,
                     FieldList<Dimension, typename Dimension::FifthRankTensor>& hessD) {
  // Make sure nodelists are of the correct size
  const auto numNodeLists = A.size();
  REQUIRE(volume.size() == numNodeLists);
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  switch (correctionOrder) {
  case CRKOrder::CubicOrder:
    REQUIRE(D.size() == numNodeLists);
    REQUIRE(gradD.size() == numNodeLists);
    REQUIRE(hessD.size() == numNodeLists);
  case CRKOrder::QuadraticOrder:
    REQUIRE(C.size() == numNodeLists);
    REQUIRE(gradC.size() == numNodeLists);
    REQUIRE(hessC.size() == numNodeLists);
  case CRKOrder::LinearOrder:
    REQUIRE(B.size() == numNodeLists);
    REQUIRE(gradB.size() == numNodeLists);
    REQUIRE(hessB.size() == numNodeLists);
  case CRKOrder::ZerothOrder:
    REQUIRE(A.size() == numNodeLists);
    REQUIRE(gradA.size() == numNodeLists);
    REQUIRE(hessA.size() == numNodeLists);
  }
  
  // Compute things point by point
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = connectivityMap.numNodes(nodeListi);
    
    for (auto nodei = 0; nodei < numNodes; ++nodei) {
      // Get list of moments
      RKMomentValues<Dimension> moments;
      
      // Calculate M values
      const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, nodei);
      for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
        for (auto nodej : connectivity[nodeListj]) {
          // Get kernel values
          const auto eta = ;
          const auto gradEta = ;
          const auto hessEta = ;
          const auto Wij = ;
          const auto gradWij = ;
          const auto hessWij = ;
          const auto vj = volume(nodeListj, nodej);

          switch (correctionOrder) {
          case CRKOrder::ZerothOrder:
            addToMomentsZeroth(eta, gradEta, hessEta, Wij, gradWij, hessWij, vj, moments);
            break;
          case CRKOrder::LinearOrder:
            addToMomentsLinear(eta, gradEta, hessEta, Wij, gradWij, hessWij, vj, moments);
            break;
          case CRKOrder::QuadraticOrder:
            addToMomentsQuadratic(eta, gradEta, hessEta, Wij, gradWij, hessWij, vj, moments);
            break;
          case CRKOrder::CubicOrder:
            addToMomentsCubic(eta, gradEta, hessEta, Wij, gradWij, hessWij, vj, moments);
            break;
          }
        }
      }

      // Calculate conditions
      switch (correctionOrder) {
      case CRKOrder::ZerothOrder:
        computeZerothCorrections(moments,
                                 A(nodeListi, nodei), gradA(nodeListi, nodei), hessA(nodeListi, nodei));
        break;
      case CRKOrder::LinearOrder:
        computeLinearCorrections(moments,
                                 A(nodeListi, nodei), B(nodeListi, nodei),
                                 gradA(nodeListi, nodei), gradB(nodeListi, nodei),
                                 hessA(nodeListi, nodei), hessB(nodeListi, nodei));
          break;
      case CRKOrder::QuadraticOrder:
        computeQuadraticCorrections(moments,
                                    A(nodeListi, nodei), B(nodeListi, nodei), C(nodeListi, nodei),
                                    gradA(nodeListi, nodei), gradB(nodeListi, nodei), gradC(nodeListi, nodei),
                                    hessA(nodeListi, nodei), hessB(nodeListi, nodei), hessC(nodeListi, nodei));
                                    
        break;
      case CRKOrder::CubicOrder:
        computeCubicCorrections(moments,
                                A(nodeListi, nodei), B(nodeListi, nodei), C(nodeListi, nodei), D(nodeListi, nodei),
                                gradA(nodeListi, nodei), gradB(nodeListi, nodei), gradC(nodeListi, nodei), gradD(nodeListi, nodei),
                                hessA(nodeListi, nodei), hessB(nodeListi, nodei), hessC(nodeListi, nodei), hessD(nodeListi, nodei));
        break;
      }
    } // nodei
  } // nodeListi
} // computeRKCorrections

} // end namespace spheral
