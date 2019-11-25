//------------------------------------------------------------------------------
// Compute the RK correction terms
//------------------------------------------------------------------------------
#include "computeRKCorrections.hh"

#include <algorithm>
#include <array>
#include "Eigen/Dense"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

namespace {

//------------------------------------------------------------------------------
// Helper class for holding the moments
//------------------------------------------------------------------------------
template<typename Dimension>
struct RKMomentValues {
  typedef double Scalar;
  typedef std::array<double, Dimension::nDim> Vector;
  typedef std::array<double, Dimension::nDim*Dimension::nDim> Tensor;
  typedef std::array<double, Dimension::nDim*Dimension::nDim*Dimension::nDim> Tensor3;
  typedef std::array<double, Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim> Tensor4;
  typedef std::array<double, Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim> Tensor5;
  typedef std::array<double, Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim> Tensor6;
  typedef std::array<double, Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim> Tensor7;
  typedef std::array<double, Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim*Dimension::nDim> Tensor8;

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

  // Zero out all values
  inline void zeroZeroth() {
    m0 = 0.0;
    std::fill(std::begin(dm0), std::end(dm0), 0.0);
    std::fill(std::begin(ddm0), std::end(ddm0), 0.0);
  }
  inline void zeroLinear() {
    zeroZeroth();
    std::fill(std::begin(m1), std::end(m1), 0.0);
    std::fill(std::begin(m2), std::end(m2), 0.0);
    std::fill(std::begin(dm1), std::end(dm1), 0.0);
    std::fill(std::begin(dm2), std::end(dm2), 0.0);
    std::fill(std::begin(ddm1), std::end(ddm1), 0.0);
    std::fill(std::begin(ddm2), std::end(ddm2), 0.0);
  }
  inline void zeroQuadratic() {
    zeroLinear();
    std::fill(std::begin(m3), std::end(m3), 0.0);
    std::fill(std::begin(m4), std::end(m4), 0.0);
    std::fill(std::begin(dm3), std::end(dm3), 0.0);
    std::fill(std::begin(dm4), std::end(dm4), 0.0);
    std::fill(std::begin(ddm3), std::end(ddm3), 0.0);
    std::fill(std::begin(ddm4), std::end(ddm4), 0.0);
  }
  inline void zeroCubic() {
    zeroQuadratic();
    std::fill(std::begin(m5), std::end(m5), 0.0);
    std::fill(std::begin(m6), std::end(m6), 0.0);
    std::fill(std::begin(dm5), std::end(dm5), 0.0);
    std::fill(std::begin(dm6), std::end(dm6), 0.0);
    std::fill(std::begin(ddm5), std::end(ddm5), 0.0);
    std::fill(std::begin(ddm6), std::end(ddm6), 0.0);
  }
  template<CRKOrder correctionOrder>
  inline
  void zero();
};

template<>
template<>
inline
void
RKMomentValues<Dim<1>>::
zero<CRKOrder::ZerothOrder>() {
  zeroZeroth();
}

template<>
template<>
inline
void
RKMomentValues<Dim<1>>::
zero<CRKOrder::LinearOrder>() {
  zeroLinear();
}

template<>
template<>
inline
void
RKMomentValues<Dim<1>>::
zero<CRKOrder::QuadraticOrder>() {
  zeroQuadratic();
}

template<>
template<>
inline
void
RKMomentValues<Dim<1>>::
zero<CRKOrder::CubicOrder>() {
  zeroCubic();
}

template<>
template<>
inline
void
RKMomentValues<Dim<2>>::
zero<CRKOrder::ZerothOrder>() {
  zeroZeroth();
}

template<>
template<>
inline
void
RKMomentValues<Dim<2>>::
zero<CRKOrder::LinearOrder>() {
  zeroLinear();
}

template<>
template<>
inline
void
RKMomentValues<Dim<2>>::
zero<CRKOrder::QuadraticOrder>() {
  zeroQuadratic();
}

template<>
template<>
inline
void
RKMomentValues<Dim<2>>::
zero<CRKOrder::CubicOrder>() {
  zeroCubic();
}

template<>
template<>
inline
void
RKMomentValues<Dim<3>>::
zero<CRKOrder::ZerothOrder>() {
  zeroZeroth();
}

template<>
template<>
inline
void
RKMomentValues<Dim<3>>::
zero<CRKOrder::LinearOrder>() {
  zeroLinear();
}

template<>
template<>
inline
void
RKMomentValues<Dim<3>>::
zero<CRKOrder::QuadraticOrder>() {
  zeroQuadratic();
}

template<>
template<>
inline
void
RKMomentValues<Dim<3>>::
zero<CRKOrder::CubicOrder>() {
  zeroCubic();
}

//------------------------------------------------------------------------------
// Add value to moments
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
inline
void
addToMoments(const typename Dimension::Vector& g,
             const typename Dimension::Tensor& dg,
             const typename Dimension::Scalar& w,
             const typename Dimension::Vector& dw,
             const typename Dimension::Tensor& ddw,
             const typename Dimension::Scalar& v,
             RKMomentValues<Dimension>& mom);

// d=1, o=0
template<>
inline
void
addToMoments<Dim<1>, CRKOrder::ZerothOrder>(const Dim<1>::Vector& g,
                                            const Dim<1>::Tensor& dg,
                                            const Dim<1>::Scalar& w,
                                            const Dim<1>::Vector& dw,
                                            const Dim<1>::Tensor& ddw,
                                            const Dim<1>::Scalar& v,
                                            RKMomentValues<Dim<1>>& mom) {
  const auto dim = Dim<1>::nDim;
  auto& m0 = mom.m0;
  auto& dm0 = mom.dm0;
  auto& ddm0 = mom.ddm0;
  const auto k1 = 0;
  const auto k2 = 0;
  
  // Moments
  m0 += v*w;
  // Gradients
  dm0[k1] += v*dw[k1];
  // Hessians
  ddm0[k1 + k2] += v*ddw[k1 + k2];
}
// d=2, o=0
template<>
inline
void
addToMoments<Dim<2>, CRKOrder::ZerothOrder>(const Dim<2>::Vector& g,
                                            const Dim<2>::Tensor& dg,
                                            const Dim<2>::Scalar& w,
                                            const Dim<2>::Vector& dw,
                                            const Dim<2>::Tensor& ddw,
                                            const Dim<2>::Scalar& v,
                                            RKMomentValues<Dim<2>>& mom) {
  const auto dim = Dim<2>::nDim;
  auto& m0 = mom.m0;
  auto& dm0 = mom.dm0;
  auto& ddm0 = mom.ddm0;
  
  // Moments
  m0 += v*w;
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    dm0[k1] += v*dw[k1];
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      ddm0[2*k1 + k2] += v*ddw[2*k1 + k2];
    }
  }
}
// d=3, o=0
template<>
inline
void
addToMoments<Dim<3>, CRKOrder::ZerothOrder>(const Dim<3>::Vector& g,
                                            const Dim<3>::Tensor& dg,
                                            const Dim<3>::Scalar& w,
                                            const Dim<3>::Vector& dw,
                                            const Dim<3>::Tensor& ddw,
                                            const Dim<3>::Scalar& v,
                                            RKMomentValues<Dim<3>>& mom) {
  const auto dim = Dim<3>::nDim;
  auto& m0 = mom.m0;
  auto& dm0 = mom.dm0;
  auto& ddm0 = mom.ddm0;
  
  // Moments
  m0 += v*w;
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    dm0[k1] += v*dw[k1];
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      ddm0[3*k1 + k2] += v*ddw[3*k1 + k2];
    }
  }
}
// d=1, o=1
template<>
inline
void
addToMoments<Dim<1>, CRKOrder::LinearOrder>(const Dim<1>::Vector& g,
                                            const Dim<1>::Tensor& dg,
                                            const Dim<1>::Scalar& w,
                                            const Dim<1>::Vector& dw,
                                            const Dim<1>::Tensor& ddw,
                                            const Dim<1>::Scalar& v,
                                            RKMomentValues<Dim<1>>& mom) {
  const auto dim = Dim<1>::nDim;
  auto& m1 = mom.m1;
  auto& dm1 = mom.dm1;
  auto& ddm1 = mom.ddm1;
  auto& m2 = mom.m2;
  auto& dm2 = mom.dm2;
  auto& ddm2 = mom.ddm2;
  const auto k1 = 0;
  const auto k2 = 0;
  const auto q1 = 0;
  const auto q2 = 0;

  // Previous moments
  addToMoments<Dim<1>, CRKOrder::ZerothOrder>(g, dg, 
                                              w, dw, ddw, v,
                                              mom);
  
  // Moments
  m1[q1] += v*w*g[q1];
  m2[q1 + q2] += v*w*g[q1]*g[q2];

  // Gradients
  dm1[k1 + q1] += v*dw[k1]*g[q1] + v*w*dg[k1 + q1];
  dm2[k1 + q1 + q2] += v*dw[k1]*g[q1]*g[q2] + v*w*g[q2]*dg[k1 + q1] + v*w*g[q1]*dg[k1 + q2];
  
  // Hessians
  ddm1[k1 + k2 + q1] += v*ddw[k1 + k2]*g[q1] + v*dw[k2]*dg[k1 + q1] + v*dw[k1]*dg[k2 + q1];
  ddm2[k1 + k2 + q1 + q2] += v*ddw[k1 + k2]*g[q1]*g[q2] + v*dw[k2]*g[q2]*dg[k1 + q1] + v*dw[k1]*g[q2]*dg[k2 + q1] + v*dw[k2]*g[q1]*dg[k1 + q2] + v*w*dg[k2 + q1]*dg[k1 + q2] + v*dw[k1]*g[q1]*dg[k2 + q2] + v*w*dg[k1 + q1]*dg[k2 + q2];
}
// d=2, o=1
template<>
inline
void
addToMoments<Dim<2>, CRKOrder::LinearOrder>(const Dim<2>::Vector& g,
                                            const Dim<2>::Tensor& dg,
                                            const Dim<2>::Scalar& w,
                                            const Dim<2>::Vector& dw,
                                            const Dim<2>::Tensor& ddw,
                                            const Dim<2>::Scalar& v,
                                            RKMomentValues<Dim<2>>& mom) {
  const auto dim = Dim<2>::nDim;
  auto& m1 = mom.m1;
  auto& dm1 = mom.dm1;
  auto& ddm1 = mom.ddm1;
  auto& m2 = mom.m2;
  auto& dm2 = mom.dm2;
  auto& ddm2 = mom.ddm2;
  
  // Previous moments
  addToMoments<Dim<2>, CRKOrder::ZerothOrder>(g, dg,
                                              w, dw, ddw, v,
                                              mom);
  
  // Moments
  for (auto q1 = 0; q1 < dim; ++q1) {
    m1[q1] += v*w*g[q1];
    for (auto q2 = 0; q2 < dim; ++q2) {
      m2[2*q1 + q2] += v*w*g[q1]*g[q2];
    }
  }
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto q1 = 0; q1 < dim; ++q1) {
      dm1[2*q1 + k1] += v*dw[k1]*g[q1] + v*w*dg[k1 + 2*q1];
      for (auto q2 = q1; q2 < dim; ++q2) {
        dm2[4*q1 + 2*q2 + k1] += v*dw[k1]*g[q1]*g[q2] + v*w*g[q2]*dg[k1 + 2*q1] + v*w*g[q1]*dg[k1 + 2*q2];
      }
    }
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        ddm1[4*q1 + 2*k1 + k2] += v*ddw[2*k1 + k2]*g[q1] + v*dw[k2]*dg[k1 + 2*q1] + v*dw[k1]*dg[k2 + 2*q1];
        for (auto q2 = q1; q2 < dim; ++q2) {
          ddm2[8*q1 + 4*q2 + 2*k1 + k2] += v*ddw[2*k1 + k2]*g[q1]*g[q2] + v*dw[k2]*g[q2]*dg[k1 + 2*q1] + v*dw[k1]*g[q2]*dg[k2 + 2*q1] + v*dw[k2]*g[q1]*dg[k1 + 2*q2] + v*w*dg[k2 + 2*q1]*dg[k1 + 2*q2] + v*dw[k1]*g[q1]*dg[k2 + 2*q2] + v*w*dg[k1 + 2*q1]*dg[k2 + 2*q2];
        }
      }
    }
  }
}
// d=3, o=1
template<>
inline
void
addToMoments<Dim<3>, CRKOrder::LinearOrder>(const Dim<3>::Vector& g,
                                            const Dim<3>::Tensor& dg,
                                            const Dim<3>::Scalar& w,
                                            const Dim<3>::Vector& dw,
                                            const Dim<3>::Tensor& ddw,
                                            const Dim<3>::Scalar& v,
                                            RKMomentValues<Dim<3>>& mom) {
  const auto dim = Dim<3>::nDim;
  auto& m1 = mom.m1;
  auto& dm1 = mom.dm1;
  auto& ddm1 = mom.ddm1;
  auto& m2 = mom.m2;
  auto& dm2 = mom.dm2;
  auto& ddm2 = mom.ddm2;
  
  // Previous moments
  addToMoments<Dim<3>, CRKOrder::ZerothOrder>(g, dg, 
                                              w, dw, ddw, v,
                                              mom);
  
  // Moments
  for (auto q1 = 0; q1 < dim; ++q1) {
    m1[q1] += v*w*g[q1];
    for (auto q2 = q1; q2 < dim; ++q2) {
      m2[3*q1 + q2] += v*w*g[q1]*g[q2];
    }
  }
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto q1 = 0; q1 < dim; ++q1) {
      dm1[3*q1 + k1] += v*dw[k1]*g[q1] + v*w*dg[k1 + 3*q1];
      for (auto q2 = q1; q2 < dim; ++q2) {
        dm2[9*q1 + 3*q2 + k1] += v*dw[k1]*g[q1]*g[q2] + v*w*g[q2]*dg[k1 + 3*q1] + v*w*g[q1]*dg[k1 + 3*q2];
      }
    }
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        ddm1[9*q1 + 3*k1 + k2] += v*ddw[3*k1 + k2]*g[q1] + v*dw[k2]*dg[k1 + 3*q1] + v*dw[k1]*dg[k2 + 3*q1];
        for (auto q2 = q1; q2 < dim; ++q2) {
          ddm2[27*q1 + 9*q2 + 3*k1 + k2] += v*ddw[3*k1 + k2]*g[q1]*g[q2] + v*dw[k2]*g[q2]*dg[k1 + 3*q1] + v*dw[k1]*g[q2]*dg[k2 + 3*q1] + v*dw[k2]*g[q1]*dg[k1 + 3*q2] + v*w*dg[k2 + 3*q1]*dg[k1 + 3*q2] + v*dw[k1]*g[q1]*dg[k2 + 3*q2] + v*w*dg[k1 + 3*q1]*dg[k2 + 3*q2];
        }
      }
    }
  }
}
// d=1, o=2
template<>
inline
void
addToMoments<Dim<1>, CRKOrder::QuadraticOrder>(const Dim<1>::Vector& g,
                                               const Dim<1>::Tensor& dg,
                                               const Dim<1>::Scalar& w,
                                               const Dim<1>::Vector& dw,
                                               const Dim<1>::Tensor& ddw,
                                               const Dim<1>::Scalar& v,
                                               RKMomentValues<Dim<1>>& mom) {
  const auto dim = Dim<1>::nDim;
  auto& m3 = mom.m3;
  auto& m4 = mom.m4;
  auto& dm3 = mom.dm3;
  auto& dm4 = mom.dm4;
  auto& ddm3 = mom.ddm3;
  auto& ddm4 = mom.ddm4;
  const auto k1 = 0;
  const auto k2 = 0;
  const auto q1 = 0;
  const auto q2 = 0;
  const auto q3 = 0;
  const auto q4 = 0;

  // Previous moments
  addToMoments<Dim<1>, CRKOrder::LinearOrder>(g, dg, 
                                              w, dw, ddw, v,
                                              mom);
  
  // Moments
  m3[q1 + q2 + q3] += v*w*g[q1]*g[q2]*g[q3];
  m4[q1 + q2 + q3 + q4] += v*w*g[q1]*g[q2]*g[q3]*g[q4];

  // Gradients
  dm3[k1 + q1 + q2 + q3] += v*dw[k1]*g[q1]*g[q2]*g[q3] + v*w*g[q2]*g[q3]*dg[k1 + q1] + v*w*g[q1]*g[q3]*dg[k1 + q2] + v*w*g[q1]*g[q2]*dg[k1 + q3];
  dm4[k1 + q1 + q2 + q3 + q4] += v*dw[k1]*g[q1]*g[q2]*g[q3]*g[q4] + v*w*g[q2]*g[q3]*g[q4]*dg[k1 + q1] + v*w*g[q1]*g[q3]*g[q4]*dg[k1 + q2] + v*w*g[q1]*g[q2]*g[q4]*dg[k1 + q3] + v*w*g[q1]*g[q2]*g[q3]*dg[k1 + q4];
  
  // Hessians
  ddm3[k1 + k2 + q1 + q2 + q3] += v*ddw[k1 + k2]*g[q1]*g[q2]*g[q3] + v*dw[k2]*g[q2]*g[q3]*dg[k1 + q1] + v*dw[k1]*g[q2]*g[q3]*dg[k2 + q1] + v*dw[k2]*g[q1]*g[q3]*dg[k1 + q2] + v*w*g[q3]*dg[k2 + q1]*dg[k1 + q2] + v*dw[k1]*g[q1]*g[q3]*dg[k2 + q2] + v*w*g[q3]*dg[k1 + q1]*dg[k2 + q2] + v*dw[k2]*g[q1]*g[q2]*dg[k1 + q3] + v*w*g[q2]*dg[k2 + q1]*dg[k1 + q3] + v*w*g[q1]*dg[k2 + q2]*dg[k1 + q3] + v*dw[k1]*g[q1]*g[q2]*dg[k2 + q3] + v*w*g[q2]*dg[k1 + q1]*dg[k2 + q3] + v*w*g[q1]*dg[k1 + q2]*dg[k2 + q3];
  ddm4[k1 + k2 + q1 + q2 + q3 + q4] += v*ddw[k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4] + v*dw[k2]*g[q2]*g[q3]*g[q4]*dg[k1 + q1] + v*dw[k1]*g[q2]*g[q3]*g[q4]*dg[k2 + q1] + v*dw[k2]*g[q1]*g[q3]*g[q4]*dg[k1 + q2] + v*w*g[q3]*g[q4]*dg[k2 + q1]*dg[k1 + q2] + v*dw[k1]*g[q1]*g[q3]*g[q4]*dg[k2 + q2] + v*w*g[q3]*g[q4]*dg[k1 + q1]*dg[k2 + q2] + v*dw[k2]*g[q1]*g[q2]*g[q4]*dg[k1 + q3] + v*w*g[q2]*g[q4]*dg[k2 + q1]*dg[k1 + q3] + v*w*g[q1]*g[q4]*dg[k2 + q2]*dg[k1 + q3] + v*dw[k1]*g[q1]*g[q2]*g[q4]*dg[k2 + q3] + v*w*g[q2]*g[q4]*dg[k1 + q1]*dg[k2 + q3] + v*w*g[q1]*g[q4]*dg[k1 + q2]*dg[k2 + q3] + v*dw[k2]*g[q1]*g[q2]*g[q3]*dg[k1 + q4] + v*w*g[q2]*g[q3]*dg[k2 + q1]*dg[k1 + q4] + v*w*g[q1]*g[q3]*dg[k2 + q2]*dg[k1 + q4] + v*w*g[q1]*g[q2]*dg[k2 + q3]*dg[k1 + q4] + v*dw[k1]*g[q1]*g[q2]*g[q3]*dg[k2 + q4] + v*w*g[q2]*g[q3]*dg[k1 + q1]*dg[k2 + q4] + v*w*g[q1]*g[q3]*dg[k1 + q2]*dg[k2 + q4] + v*w*g[q1]*g[q2]*dg[k1 + q3]*dg[k2 + q4];
}
// d=2, o=2
template<>
inline
void
addToMoments<Dim<2>, CRKOrder::QuadraticOrder>(const Dim<2>::Vector& g,
                                               const Dim<2>::Tensor& dg,
                                               const Dim<2>::Scalar& w,
                                               const Dim<2>::Vector& dw,
                                               const Dim<2>::Tensor& ddw,
                                               const Dim<2>::Scalar& v,
                                               RKMomentValues<Dim<2>>& mom) {
  const auto dim = Dim<2>::nDim;
  auto& m3 = mom.m3;
  auto& m4 = mom.m4;
  auto& dm3 = mom.dm3;
  auto& dm4 = mom.dm4;
  auto& ddm3 = mom.ddm3;
  auto& ddm4 = mom.ddm4;
  
  // Previous moments
  addToMoments<Dim<2>, CRKOrder::LinearOrder>(g, dg, 
                                              w, dw, ddw, v,
                                              mom);
  
  // Moments
  for (auto q1 = 0; q1 < dim; ++q1) {
    for (auto q2 = q1; q2 < dim; ++q2) {
      for (auto q3 = q2; q3 < dim; ++q3) {
        m3[4*q1 + 2*q2 + q3] += v*w*g[q1]*g[q2]*g[q3];
        for (auto q4 = q3; q4 < dim; ++q4) {
          m4[8*q1 + 4*q2 + 2*q3 + q4] += v*w*g[q1]*g[q2]*g[q3]*g[q4];
        }
      }
    }
  }
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto q1 = 0; q1 < dim; ++q1) {
      for (auto q2 = q1; q2 < dim; ++q2) {
        for (auto q3 = q2; q3 < dim; ++q3) {
          dm3[8*q1 + 4*q2 + 2*q3 + k1] += v*dw[k1]*g[q1]*g[q2]*g[q3] + v*w*g[q2]*g[q3]*dg[k1 + 2*q1] + v*w*g[q1]*g[q3]*dg[k1 + 2*q2] + v*w*g[q1]*g[q2]*dg[k1 + 2*q3];
          for (auto q4 = q3; q4 < dim; ++q4) {
            dm4[16*q1 + 8*q2 + 4*q3 + 2*q4 + k1] += v*dw[k1]*g[q1]*g[q2]*g[q3]*g[q4] + v*w*g[q2]*g[q3]*g[q4]*dg[k1 + 2*q1] + v*w*g[q1]*g[q3]*g[q4]*dg[k1 + 2*q2] + v*w*g[q1]*g[q2]*g[q4]*dg[k1 + 2*q3] + v*w*g[q1]*g[q2]*g[q3]*dg[k1 + 2*q4];
          }
        }
      }
    }
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        for (auto q2 = q1; q2 < dim; ++q2) {
          for (auto q3 = q2; q3 < dim; ++q3) {
            ddm3[16*q1 + 8*q2 + 4*q3 + 2*k1 + k2] += v*ddw[2*k1 + k2]*g[q1]*g[q2]*g[q3] + v*dw[k2]*g[q2]*g[q3]*dg[k1 + 2*q1] + v*dw[k1]*g[q2]*g[q3]*dg[k2 + 2*q1] + v*dw[k2]*g[q1]*g[q3]*dg[k1 + 2*q2] + v*w*g[q3]*dg[k2 + 2*q1]*dg[k1 + 2*q2] + v*dw[k1]*g[q1]*g[q3]*dg[k2 + 2*q2] + v*w*g[q3]*dg[k1 + 2*q1]*dg[k2 + 2*q2] + v*dw[k2]*g[q1]*g[q2]*dg[k1 + 2*q3] + v*w*g[q2]*dg[k2 + 2*q1]*dg[k1 + 2*q3] + v*w*g[q1]*dg[k2 + 2*q2]*dg[k1 + 2*q3] + v*dw[k1]*g[q1]*g[q2]*dg[k2 + 2*q3] + v*w*g[q2]*dg[k1 + 2*q1]*dg[k2 + 2*q3] + v*w*g[q1]*dg[k1 + 2*q2]*dg[k2 + 2*q3];
            for (auto q4 = q3; q4 < dim; ++q4) {
              ddm4[32*q1 +16*q2 + 8*q3 + 4*q4 + 2*k1] += v*ddw[2*k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4] + v*dw[k2]*g[q2]*g[q3]*g[q4]*dg[k1 + 2*q1] + v*dw[k1]*g[q2]*g[q3]*g[q4]*dg[k2 + 2*q1] + v*dw[k2]*g[q1]*g[q3]*g[q4]*dg[k1 + 2*q2] + v*w*g[q3]*g[q4]*dg[k2 + 2*q1]*dg[k1 + 2*q2] + v*dw[k1]*g[q1]*g[q3]*g[q4]*dg[k2 + 2*q2] + v*w*g[q3]*g[q4]*dg[k1 + 2*q1]*dg[k2 + 2*q2] + v*dw[k2]*g[q1]*g[q2]*g[q4]*dg[k1 + 2*q3] + v*w*g[q2]*g[q4]*dg[k2 + 2*q1]*dg[k1 + 2*q3] + v*w*g[q1]*g[q4]*dg[k2 + 2*q2]*dg[k1 + 2*q3] + v*dw[k1]*g[q1]*g[q2]*g[q4]*dg[k2 + 2*q3] + v*w*g[q2]*g[q4]*dg[k1 + 2*q1]*dg[k2 + 2*q3] + v*w*g[q1]*g[q4]*dg[k1 + 2*q2]*dg[k2 + 2*q3] + v*dw[k2]*g[q1]*g[q2]*g[q3]*dg[k1 + 2*q4] + v*w*g[q2]*g[q3]*dg[k2 + 2*q1]*dg[k1 + 2*q4] + v*w*g[q1]*g[q3]*dg[k2 + 2*q2]*dg[k1 + 2*q4] + v*w*g[q1]*g[q2]*dg[k2 + 2*q3]*dg[k1 + 2*q4] + v*dw[k1]*g[q1]*g[q2]*g[q3]*dg[k2 + 2*q4] + v*w*g[q2]*g[q3]*dg[k1 + 2*q1]*dg[k2 + 2*q4] + v*w*g[q1]*g[q3]*dg[k1 + 2*q2]*dg[k2 + 2*q4] + v*w*g[q1]*g[q2]*dg[k1 + 2*q3]*dg[k2 + 2*q4];
            }
          }
        }
      }
    }
  }
}
// d=3, o=2
template<>
inline
void
addToMoments<Dim<3>, CRKOrder::QuadraticOrder>(const Dim<3>::Vector& g,
                                               const Dim<3>::Tensor& dg,
                                               const Dim<3>::Scalar& w,
                                               const Dim<3>::Vector& dw,
                                               const Dim<3>::Tensor& ddw,
                                               const Dim<3>::Scalar& v,
                                               RKMomentValues<Dim<3>>& mom) {
  const auto dim = Dim<3>::nDim;
  auto& m3 = mom.m3;
  auto& m4 = mom.m4;
  auto& dm3 = mom.dm3;
  auto& dm4 = mom.dm4;
  auto& ddm3 = mom.ddm3;
  auto& ddm4 = mom.ddm4;
  
  // Previous moments
  addToMoments<Dim<3>, CRKOrder::LinearOrder>(g, dg, 
                                              w, dw, ddw, v,
                                              mom);
  
  // Moments
  for (auto q1 = 0; q1 < dim; ++q1) {
    for (auto q2 = q1; q2 < dim; ++q2) {
      for (auto q3 = q2; q3 < dim; ++q3) {
        m3[9*q1 + 3*q2 + q3] += v*w*g[q1]*g[q2]*g[q3];
        for (auto q4 = q3; q4 < dim; ++q4) {
          m4[27*q1 + 9*q2 + 3*q3 + q4] += v*w*g[q1]*g[q2]*g[q3]*g[q4];
        }
      }
    }
  }
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto q1 = 0; q1 < dim; ++q1) {
      for (auto q2 = q1; q2 < dim; ++q2) {
        for (auto q3 = q2; q3 < dim; ++q3) {
          dm3[27*q1 + 9*q2 + 3*q3 + k1] += v*dw[k1]*g[q1]*g[q2]*g[q3] + v*w*g[q2]*g[q3]*dg[k1 + 3*q1] + v*w*g[q1]*g[q3]*dg[k1 + 3*q2] + v*w*g[q1]*g[q2]*dg[k1 + 3*q3];
          for (auto q4 = q3; q4 < dim; ++q4) {
            dm4[81*q1 + 27*q2 + 9*q3 + 3*q4 + k1] += v*dw[k1]*g[q1]*g[q2]*g[q3]*g[q4] + v*w*g[q2]*g[q3]*g[q4]*dg[k1 + 3*q1] + v*w*g[q1]*g[q3]*g[q4]*dg[k1 + 3*q2] + v*w*g[q1]*g[q2]*g[q4]*dg[k1 + 3*q3] + v*w*g[q1]*g[q2]*g[q3]*dg[k1 + 3*q4];
          }
        }
      }
    }
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        for (auto q2 = q1; q2 < dim; ++q2) {
          for (auto q3 = q2; q3 < dim; ++q3) {
            ddm3[81*q1 + 27*q2 + 9*q3 + 3*k1 + k2] += v*ddw[3*k1 + k2]*g[q1]*g[q2]*g[q3] + v*dw[k2]*g[q2]*g[q3]*dg[k1 + 3*q1] + v*dw[k1]*g[q2]*g[q3]*dg[k2 + 3*q1] + v*dw[k2]*g[q1]*g[q3]*dg[k1 + 3*q2] + v*w*g[q3]*dg[k2 + 3*q1]*dg[k1 + 3*q2] + v*dw[k1]*g[q1]*g[q3]*dg[k2 + 3*q2] + v*w*g[q3]*dg[k1 + 3*q1]*dg[k2 + 3*q2] + v*dw[k2]*g[q1]*g[q2]*dg[k1 + 3*q3] + v*w*g[q2]*dg[k2 + 3*q1]*dg[k1 + 3*q3] + v*w*g[q1]*dg[k2 + 3*q2]*dg[k1 + 3*q3] + v*dw[k1]*g[q1]*g[q2]*dg[k2 + 3*q3] + v*w*g[q2]*dg[k1 + 3*q1]*dg[k2 + 3*q3] + v*w*g[q1]*dg[k1 + 3*q2]*dg[k2 + 3*q3];
            for (auto q4 = q3; q4 < dim; ++q4) {
              ddm4[243*q1 +81*q2 + 27*q3 + 9*q4 + 3*k1] += v*ddw[3*k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4] + v*dw[k2]*g[q2]*g[q3]*g[q4]*dg[k1 + 3*q1] + v*dw[k1]*g[q2]*g[q3]*g[q4]*dg[k2 + 3*q1] + v*dw[k2]*g[q1]*g[q3]*g[q4]*dg[k1 + 3*q2] + v*w*g[q3]*g[q4]*dg[k2 + 3*q1]*dg[k1 + 3*q2] + v*dw[k1]*g[q1]*g[q3]*g[q4]*dg[k2 + 3*q2] + v*w*g[q3]*g[q4]*dg[k1 + 3*q1]*dg[k2 + 3*q2] + v*dw[k2]*g[q1]*g[q2]*g[q4]*dg[k1 + 3*q3] + v*w*g[q2]*g[q4]*dg[k2 + 3*q1]*dg[k1 + 3*q3] + v*w*g[q1]*g[q4]*dg[k2 + 3*q2]*dg[k1 + 3*q3] + v*dw[k1]*g[q1]*g[q2]*g[q4]*dg[k2 + 3*q3] + v*w*g[q2]*g[q4]*dg[k1 + 3*q1]*dg[k2 + 3*q3] + v*w*g[q1]*g[q4]*dg[k1 + 3*q2]*dg[k2 + 3*q3] + v*dw[k2]*g[q1]*g[q2]*g[q3]*dg[k1 + 3*q4] + v*w*g[q2]*g[q3]*dg[k2 + 3*q1]*dg[k1 + 3*q4] + v*w*g[q1]*g[q3]*dg[k2 + 3*q2]*dg[k1 + 3*q4] + v*w*g[q1]*g[q2]*dg[k2 + 3*q3]*dg[k1 + 3*q4] + v*dw[k1]*g[q1]*g[q2]*g[q3]*dg[k2 + 3*q4] + v*w*g[q2]*g[q3]*dg[k1 + 3*q1]*dg[k2 + 3*q4] + v*w*g[q1]*g[q3]*dg[k1 + 3*q2]*dg[k2 + 3*q4] + v*w*g[q1]*g[q2]*dg[k1 + 3*q3]*dg[k2 + 3*q4];
            }
          }
        }
      }
    }
  }
}
// d=1, o=3
template<>
inline
void
addToMoments<Dim<1>, CRKOrder::CubicOrder>(const Dim<1>::Vector& g,
                                           const Dim<1>::Tensor& dg,
                                           const Dim<1>::Scalar& w,
                                           const Dim<1>::Vector& dw,
                                           const Dim<1>::Tensor& ddw,
                                           const Dim<1>::Scalar& v,
                                           RKMomentValues<Dim<1>>& mom) {
  const auto dim = Dim<1>::nDim;
  auto& m5 = mom.m5;
  auto& m6 = mom.m6;
  auto& dm5 = mom.dm5;
  auto& dm6 = mom.dm6;
  auto& ddm5 = mom.ddm5;
  auto& ddm6 = mom.ddm6;
  const auto k1 = 0;
  const auto k2 = 0;
  const auto q1 = 0;
  const auto q2 = 0;
  const auto q3 = 0;
  const auto q4 = 0;
  const auto q5 = 0;
  const auto q6 = 0;

  // Previous moments
  addToMoments<Dim<1>, CRKOrder::QuadraticOrder>(g, dg, 
                                                 w, dw, ddw, v,
                                                 mom);
  
  // Moments
  m5[q1 + q2 + q3 + q4 + q5] += v*w*g[q1]*g[q2]*g[q3]*g[q4]*g[q5];
  m6[q1 + q2 + q3 + q4 + q4 + q5 + q6] += v*w*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6];

  // Gradients
  dm5[k1 + q1 + q2 + q3 + q4 + q5] += v*(dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + q2] + w*g[q2]*(g[q1]*g[q4]*g[q5]*dg[k1 + q3] + g[q3]*(g[q1]*g[q5]*dg[k1 + q4] + g[q4]*(g[q5]*dg[k1 + q1] + g[q1]*dg[k1 + q5]))));
  dm6[k1 + q1 + q2 + q3 + q4 + q4 + q5 + q6] += v*(dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + q2] + w*g[q2]*(g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + q3] + g[q3]*(g[q1]*g[q5]*g[q6]*dg[k1 + q4] + g[q4]*(g[q5]*g[q6]*dg[k1 + q1] + g[q1]*g[q6]*dg[k1 + q5] + g[q1]*g[q5]*dg[k1 + q6]))));
  
  // Hessians
  ddm5[k1 + k2 + q1 + q2 + q3 + q4 + q5] += v*(ddw[k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5] + dw[k1]*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + q1] + w*g[q3]*g[q4]*g[q5]*dg[k2 + q1]*dg[k1 + q2] + dw[k1]*g[q1]*g[q3]*g[q4]*g[q5]*dg[k2 + q2] + w*g[q3]*g[q4]*g[q5]*dg[k1 + q1]*dg[k2 + q2] + w*g[q2]*g[q4]*g[q5]*dg[k2 + q1]*dg[k1 + q3] + w*g[q1]*g[q4]*g[q5]*dg[k2 + q2]*dg[k1 + q3] + dw[k1]*g[q1]*g[q2]*g[q4]*g[q5]*dg[k2 + q3] + w*g[q2]*g[q4]*g[q5]*dg[k1 + q1]*dg[k2 + q3] + w*g[q1]*g[q4]*g[q5]*dg[k1 + q2]*dg[k2 + q3] + w*g[q2]*g[q3]*g[q5]*dg[k2 + q1]*dg[k1 + q4] + w*g[q1]*g[q3]*g[q5]*dg[k2 + q2]*dg[k1 + q4] + w*g[q1]*g[q2]*g[q5]*dg[k2 + q3]*dg[k1 + q4] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q5]*dg[k2 + q4] + w*g[q2]*g[q3]*g[q5]*dg[k1 + q1]*dg[k2 + q4] + w*g[q1]*g[q3]*g[q5]*dg[k1 + q2]*dg[k2 + q4] + w*g[q1]*g[q2]*g[q5]*dg[k1 + q3]*dg[k2 + q4] + w*g[q2]*g[q3]*g[q4]*dg[k2 + q1]*dg[k1 + q5] + w*g[q1]*g[q3]*g[q4]*dg[k2 + q2]*dg[k1 + q5] + w*g[q1]*g[q2]*g[q4]*dg[k2 + q3]*dg[k1 + q5] + w*g[q1]*g[q2]*g[q3]*dg[k2 + q4]*dg[k1 + q5] + dw[k2]*(g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + q2] + g[q2]*(g[q1]*g[q4]*g[q5]*dg[k1 + q3] + g[q3]*(g[q4]*g[q5]*dg[k1 + q1] + g[q1]*g[q5]*dg[k1 + q4] + g[q1]*g[q4]*dg[k1 + q5]))) + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*dg[k2 + q5] + w*g[q2]*g[q3]*g[q4]*dg[k1 + q1]*dg[k2 + q5] + w*g[q1]*g[q3]*g[q4]*dg[k1 + q2]*dg[k2 + q5] + w*g[q1]*g[q2]*g[q4]*dg[k1 + q3]*dg[k2 + q5] + w*g[q1]*g[q2]*g[q3]*dg[k1 + q4]*dg[k2 + q5]);
  ddm6[k1 + k2 + q1 + q2 + q3 + q4 + q4 + q5 + q6] += v*(ddw[k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6] + dw[k1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + q1] + w*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + q1]*dg[k1 + q2] + dw[k1]*g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + q2] + w*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + q1]*dg[k2 + q2] + w*g[q2]*g[q4]*g[q5]*g[q6]*dg[k2 + q1]*dg[k1 + q3] + w*g[q1]*g[q4]*g[q5]*g[q6]*dg[k2 + q2]*dg[k1 + q3] + dw[k1]*g[q1]*g[q2]*g[q4]*g[q5]*g[q6]*dg[k2 + q3] + w*g[q2]*g[q4]*g[q5]*g[q6]*dg[k1 + q1]*dg[k2 + q3] + w*g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + q2]*dg[k2 + q3] + w*g[q2]*g[q3]*g[q5]*g[q6]*dg[k2 + q1]*dg[k1 + q4] + w*g[q1]*g[q3]*g[q5]*g[q6]*dg[k2 + q2]*dg[k1 + q4] + w*g[q1]*g[q2]*g[q5]*g[q6]*dg[k2 + q3]*dg[k1 + q4] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q5]*g[q6]*dg[k2 + q4] + w*g[q2]*g[q3]*g[q5]*g[q6]*dg[k1 + q1]*dg[k2 + q4] + w*g[q1]*g[q3]*g[q5]*g[q6]*dg[k1 + q2]*dg[k2 + q4] + w*g[q1]*g[q2]*g[q5]*g[q6]*dg[k1 + q3]*dg[k2 + q4] + w*g[q2]*g[q3]*g[q4]*g[q6]*dg[k2 + q1]*dg[k1 + q5] + w*g[q1]*g[q3]*g[q4]*g[q6]*dg[k2 + q2]*dg[k1 + q5] + w*g[q1]*g[q2]*g[q4]*g[q6]*dg[k2 + q3]*dg[k1 + q5] + w*g[q1]*g[q2]*g[q3]*g[q6]*dg[k2 + q4]*dg[k1 + q5] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q6]*dg[k2 + q5] + w*g[q2]*g[q3]*g[q4]*g[q6]*dg[k1 + q1]*dg[k2 + q5] + w*g[q1]*g[q3]*g[q4]*g[q6]*dg[k1 + q2]*dg[k2 + q5] + w*g[q1]*g[q2]*g[q4]*g[q6]*dg[k1 + q3]*dg[k2 + q5] + w*g[q1]*g[q2]*g[q3]*g[q6]*dg[k1 + q4]*dg[k2 + q5] + w*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + q1]*dg[k1 + q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k2 + q2]*dg[k1 + q6] + w*g[q1]*g[q2]*g[q4]*g[q5]*dg[k2 + q3]*dg[k1 + q6] + w*g[q1]*g[q2]*g[q3]*g[q5]*dg[k2 + q4]*dg[k1 + q6] + w*g[q1]*g[q2]*g[q3]*g[q4]*dg[k2 + q5]*dg[k1 + q6] + dw[k2]*(g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + q2] + g[q2]*(g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + q3] + g[q3]*(g[q1]*g[q5]*g[q6]*dg[k1 + q4] + g[q4]*(g[q5]*g[q6]*dg[k1 + q1] + g[q1]*g[q6]*dg[k1 + q5] + g[q1]*g[q5]*dg[k1 + q6])))) + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + q6] + w*g[q2]*g[q3]*g[q4]*g[q5]*dg[k1 + q1]*dg[k2 + q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + q2]*dg[k2 + q6] + w*g[q1]*g[q2]*g[q4]*g[q5]*dg[k1 + q3]*dg[k2 + q6] + w*g[q1]*g[q2]*g[q3]*g[q5]*dg[k1 + q4]*dg[k2 + q6] + w*g[q1]*g[q2]*g[q3]*g[q4]*dg[k1 + q5]*dg[k2 + q6]);
}
// d=2, o=3
template<>
inline
void
addToMoments<Dim<2>, CRKOrder::CubicOrder>(const Dim<2>::Vector& g,
                                           const Dim<2>::Tensor& dg,
                                           const Dim<2>::Scalar& w,
                                           const Dim<2>::Vector& dw,
                                           const Dim<2>::Tensor& ddw,
                                           const Dim<2>::Scalar& v,
                                           RKMomentValues<Dim<2>>& mom) {
  const auto dim = Dim<2>::nDim;
  auto& m5 = mom.m5;
  auto& m6 = mom.m6;
  auto& dm5 = mom.dm5;
  auto& dm6 = mom.dm6;
  auto& ddm5 = mom.ddm5;
  auto& ddm6 = mom.ddm6;
  
  // Previous moments
  addToMoments<Dim<2>, CRKOrder::QuadraticOrder>(g, dg, 
                                                 w, dw, ddw, v,
                                                 mom);
  
  // Moments
  for (auto q1 = 0; q1 < dim; ++q1) {
    for (auto q2 = q1; q2 < dim; ++q2) {
      for (auto q3 = q2; q3 < dim; ++q3) {
        for (auto q4 = q3; q4 < dim; ++q4) {
          for (auto q5 = q4; q5 < dim; ++q5) {
            m5[16*q1 + 8*q2 + 4*q3 + 2*q4 + q5] += v*w*g[q1]*g[q2]*g[q3]*g[q4]*g[q5];
            for (auto q6 = q5; q6 < dim; ++q6) {
              m6[32*q1 + 16*q2 + 8*q3 + 4*q4 + 2*q5 + q6] += v*w*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6];
            }
          }
        }
      }
    }
    // Gradients
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        for (auto q2 = q1; q2 < dim; ++q2) {
          for (auto q3 = q2; q3 < dim; ++q3) {
            for (auto q4 = q3; q4 < dim; ++q4) {
              for (auto q5 = q4; q5 < dim; ++q5) {
                dm5[32*q1 + 16*q2 + 8*q3 + 4*q4 + 2*q5 + k1] += v*(dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + 2*q2] + w*g[q2]*(g[q1]*g[q4]*g[q5]*dg[k1 + 2*q3] + g[q3]*(g[q1]*g[q5]*dg[k1 + 2*q4] + g[q4]*(g[q5]*dg[k1 + 2*q1] + g[q1]*dg[k1 + 2*q5]))));
                for (auto q6 = q5; q6 < dim; ++q6) {
                  dm6[64*q1 + 32*q2 + 16*q3 + 8*q4 + 4*q5 + 2*q6 + k1] += v*(dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + 2*q2] + w*g[q2]*(g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + 2*q3] + g[q3]*(g[q1]*g[q5]*g[q6]*dg[k1 + 2*q4] + g[q4]*(g[q5]*g[q6]*dg[k1 + 2*q1] + g[q1]*g[q6]*dg[k1 + 2*q5] + g[q1]*g[q5]*dg[k1 + 2*q6]))));
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
        for (auto q1 = 0; q1 < dim; ++q1) {
          for (auto q2 = q1; q2 < dim; ++q2) {
            for (auto q3 = q2; q3 < dim; ++q3) {
              for (auto q4 = q3; q4 < dim; ++q4) {
                for (auto q5 = q4; q5 < dim; ++q5) {
                  ddm5[64*q1 + 32*q2 + 16*q3 + 8*q4 + 4*q5 + 2*k1 + k2] += v*(ddw[2*k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5] + dw[k1]*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + 2*q1] + w*g[q3]*g[q4]*g[q5]*dg[k2 + 2*q1]*dg[k1 + 2*q2] + dw[k1]*g[q1]*g[q3]*g[q4]*g[q5]*dg[k2 + 2*q2] + w*g[q3]*g[q4]*g[q5]*dg[k1 + 2*q1]*dg[k2 + 2*q2] + w*g[q2]*g[q4]*g[q5]*dg[k2 + 2*q1]*dg[k1 + 2*q3] + w*g[q1]*g[q4]*g[q5]*dg[k2 + 2*q2]*dg[k1 + 2*q3] + dw[k1]*g[q1]*g[q2]*g[q4]*g[q5]*dg[k2 + 2*q3] + w*g[q2]*g[q4]*g[q5]*dg[k1 + 2*q1]*dg[k2 + 2*q3] + w*g[q1]*g[q4]*g[q5]*dg[k1 + 2*q2]*dg[k2 + 2*q3] + w*g[q2]*g[q3]*g[q5]*dg[k2 + 2*q1]*dg[k1 + 2*q4] + w*g[q1]*g[q3]*g[q5]*dg[k2 + 2*q2]*dg[k1 + 2*q4] + w*g[q1]*g[q2]*g[q5]*dg[k2 + 2*q3]*dg[k1 + 2*q4] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q5]*dg[k2 + 2*q4] + w*g[q2]*g[q3]*g[q5]*dg[k1 + 2*q1]*dg[k2 + 2*q4] + w*g[q1]*g[q3]*g[q5]*dg[k1 + 2*q2]*dg[k2 + 2*q4] + w*g[q1]*g[q2]*g[q5]*dg[k1 + 2*q3]*dg[k2 + 2*q4] + w*g[q2]*g[q3]*g[q4]*dg[k2 + 2*q1]*dg[k1 + 2*q5] + w*g[q1]*g[q3]*g[q4]*dg[k2 + 2*q2]*dg[k1 + 2*q5] + w*g[q1]*g[q2]*g[q4]*dg[k2 + 2*q3]*dg[k1 + 2*q5] + w*g[q1]*g[q2]*g[q3]*dg[k2 + 2*q4]*dg[k1 + 2*q5] + dw[k2]*(g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + 2*q2] + g[q2]*(g[q1]*g[q4]*g[q5]*dg[k1 + 2*q3] + g[q3]*(g[q4]*g[q5]*dg[k1 + 2*q1] + g[q1]*g[q5]*dg[k1 + 2*q4] + g[q1]*g[q4]*dg[k1 + 2*q5]))) + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*dg[k2 + 2*q5] + w*g[q2]*g[q3]*g[q4]*dg[k1 + 2*q1]*dg[k2 + 2*q5] + w*g[q1]*g[q3]*g[q4]*dg[k1 + 2*q2]*dg[k2 + 2*q5] + w*g[q1]*g[q2]*g[q4]*dg[k1 + 2*q3]*dg[k2 + 2*q5] + w*g[q1]*g[q2]*g[q3]*dg[k1 + 2*q4]*dg[k2 + 2*q5]);
                  for (auto q6 = q5; q6 < dim; ++q6) {
                    ddm6[128*q1 + 64*q2 + 32*q3 + 16*q4 + 8*q5 + 4*q6 + 2*k1] += v*(ddw[2*k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6] + dw[k1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + 2*q1] + w*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + 2*q1]*dg[k1 + 2*q2] + dw[k1]*g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + 2*q2] + w*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + 2*q1]*dg[k2 + 2*q2] + w*g[q2]*g[q4]*g[q5]*g[q6]*dg[k2 + 2*q1]*dg[k1 + 2*q3] + w*g[q1]*g[q4]*g[q5]*g[q6]*dg[k2 + 2*q2]*dg[k1 + 2*q3] + dw[k1]*g[q1]*g[q2]*g[q4]*g[q5]*g[q6]*dg[k2 + 2*q3] + w*g[q2]*g[q4]*g[q5]*g[q6]*dg[k1 + 2*q1]*dg[k2 + 2*q3] + w*g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + 2*q2]*dg[k2 + 2*q3] + w*g[q2]*g[q3]*g[q5]*g[q6]*dg[k2 + 2*q1]*dg[k1 + 2*q4] + w*g[q1]*g[q3]*g[q5]*g[q6]*dg[k2 + 2*q2]*dg[k1 + 2*q4] + w*g[q1]*g[q2]*g[q5]*g[q6]*dg[k2 + 2*q3]*dg[k1 + 2*q4] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q5]*g[q6]*dg[k2 + 2*q4] + w*g[q2]*g[q3]*g[q5]*g[q6]*dg[k1 + 2*q1]*dg[k2 + 2*q4] + w*g[q1]*g[q3]*g[q5]*g[q6]*dg[k1 + 2*q2]*dg[k2 + 2*q4] + w*g[q1]*g[q2]*g[q5]*g[q6]*dg[k1 + 2*q3]*dg[k2 + 2*q4] + w*g[q2]*g[q3]*g[q4]*g[q6]*dg[k2 + 2*q1]*dg[k1 + 2*q5] + w*g[q1]*g[q3]*g[q4]*g[q6]*dg[k2 + 2*q2]*dg[k1 + 2*q5] + w*g[q1]*g[q2]*g[q4]*g[q6]*dg[k2 + 2*q3]*dg[k1 + 2*q5] + w*g[q1]*g[q2]*g[q3]*g[q6]*dg[k2 + 2*q4]*dg[k1 + 2*q5] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q6]*dg[k2 + 2*q5] + w*g[q2]*g[q3]*g[q4]*g[q6]*dg[k1 + 2*q1]*dg[k2 + 2*q5] + w*g[q1]*g[q3]*g[q4]*g[q6]*dg[k1 + 2*q2]*dg[k2 + 2*q5] + w*g[q1]*g[q2]*g[q4]*g[q6]*dg[k1 + 2*q3]*dg[k2 + 2*q5] + w*g[q1]*g[q2]*g[q3]*g[q6]*dg[k1 + 2*q4]*dg[k2 + 2*q5] + w*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + 2*q1]*dg[k1 + 2*q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k2 + 2*q2]*dg[k1 + 2*q6] + w*g[q1]*g[q2]*g[q4]*g[q5]*dg[k2 + 2*q3]*dg[k1 + 2*q6] + w*g[q1]*g[q2]*g[q3]*g[q5]*dg[k2 + 2*q4]*dg[k1 + 2*q6] + w*g[q1]*g[q2]*g[q3]*g[q4]*dg[k2 + 2*q5]*dg[k1 + 2*q6] + dw[k2]*(g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + 2*q2] + g[q2]*(g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + 2*q3] + g[q3]*(g[q1]*g[q5]*g[q6]*dg[k1 + 2*q4] + g[q4]*(g[q5]*g[q6]*dg[k1 + 2*q1] + g[q1]*g[q6]*dg[k1 + 2*q5] + g[q1]*g[q5]*dg[k1 + 2*q6])))) + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + 2*q6] + w*g[q2]*g[q3]*g[q4]*g[q5]*dg[k1 + 2*q1]*dg[k2 + 2*q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + 2*q2]*dg[k2 + 2*q6] + w*g[q1]*g[q2]*g[q4]*g[q5]*dg[k1 + 2*q3]*dg[k2 + 2*q6] + w*g[q1]*g[q2]*g[q3]*g[q5]*dg[k1 + 2*q4]*dg[k2 + 2*q6] + w*g[q1]*g[q2]*g[q3]*g[q4]*dg[k1 + 2*q5]*dg[k2 + 2*q6]);
                  }
                }
              }
            }
          }
        }
      }
    }
  }
}
// d=3, o=3
template<>
inline
void
addToMoments<Dim<3>, CRKOrder::CubicOrder>(const Dim<3>::Vector& g,
                                           const Dim<3>::Tensor& dg,
                                           const Dim<3>::Scalar& w,
                                           const Dim<3>::Vector& dw,
                                           const Dim<3>::Tensor& ddw,
                                           const Dim<3>::Scalar& v,
                                           RKMomentValues<Dim<3>>& mom) {
  const auto dim = Dim<3>::nDim;
  auto& m5 = mom.m5;
  auto& m6 = mom.m6;
  auto& dm5 = mom.dm5;
  auto& dm6 = mom.dm6;
  auto& ddm5 = mom.ddm5;
  auto& ddm6 = mom.ddm6;
  
  // Previous moments
  addToMoments<Dim<3>, CRKOrder::QuadraticOrder>(g, dg, 
                                                 w, dw, ddw, v,
                                                 mom);
  
  // Moments
  for (auto q1 = 0; q1 < dim; ++q1) {
    for (auto q2 = q1; q2 < dim; ++q2) {
      for (auto q3 = q2; q3 < dim; ++q3) {
        for (auto q4 = q3; q4 < dim; ++q4) {
          for (auto q5 = q4; q5 < dim; ++q5) {
            m5[81*q1 + 27*q2 + 9*q3 + 3*q4 + q5] += v*w*g[q1]*g[q2]*g[q3]*g[q4]*g[q5];
            for (auto q6 = q5; q6 < dim; ++q6) {
              m6[243*q1 + 81*q2 + 27*q3 + 9*q4 + 3*q5 + q6] += v*w*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6];
            }
          }
        }
      }
    }
    // Gradients
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        for (auto q2 = q1; q2 < dim; ++q2) {
          for (auto q3 = q2; q3 < dim; ++q3) {
            for (auto q4 = q3; q4 < dim; ++q4) {
              for (auto q5 = q4; q5 < dim; ++q5) {
                dm5[243*q1 + 81*q2 + 27*q3 + 9*q4 + 3*q5 + k1] += v*(dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + 3*q2] + w*g[q2]*(g[q1]*g[q4]*g[q5]*dg[k1 + 3*q3] + g[q3]*(g[q1]*g[q5]*dg[k1 + 3*q4] + g[q4]*(g[q5]*dg[k1 + 3*q1] + g[q1]*dg[k1 + 3*q5]))));
                for (auto q6 = q5; q6 < dim; ++q6) {
                  dm6[729*q1 + 243*q2 + 81*q3 + 27*q4 + 9*q5 + 3*q6 + k1] += v*(dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + 3*q2] + w*g[q2]*(g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + 3*q3] + g[q3]*(g[q1]*g[q5]*g[q6]*dg[k1 + 3*q4] + g[q4]*(g[q5]*g[q6]*dg[k1 + 3*q1] + g[q1]*g[q6]*dg[k1 + 3*q5] + g[q1]*g[q5]*dg[k1 + 3*q6]))));
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
        for (auto q1 = 0; q1 < dim; ++q1) {
          for (auto q2 = q1; q2 < dim; ++q2) {
            for (auto q3 = q2; q3 < dim; ++q3) {
              for (auto q4 = q3; q4 < dim; ++q4) {
                for (auto q5 = q4; q5 < dim; ++q5) {
                  ddm5[729*q1 + 243*q2 + 81*q3 + 27*q4 + 9*q5 + 3*k1 + k2] += v*(ddw[3*k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5] + dw[k1]*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + 3*q1] + w*g[q3]*g[q4]*g[q5]*dg[k2 + 3*q1]*dg[k1 + 3*q2] + dw[k1]*g[q1]*g[q3]*g[q4]*g[q5]*dg[k2 + 3*q2] + w*g[q3]*g[q4]*g[q5]*dg[k1 + 3*q1]*dg[k2 + 3*q2] + w*g[q2]*g[q4]*g[q5]*dg[k2 + 3*q1]*dg[k1 + 3*q3] + w*g[q1]*g[q4]*g[q5]*dg[k2 + 3*q2]*dg[k1 + 3*q3] + dw[k1]*g[q1]*g[q2]*g[q4]*g[q5]*dg[k2 + 3*q3] + w*g[q2]*g[q4]*g[q5]*dg[k1 + 3*q1]*dg[k2 + 3*q3] + w*g[q1]*g[q4]*g[q5]*dg[k1 + 3*q2]*dg[k2 + 3*q3] + w*g[q2]*g[q3]*g[q5]*dg[k2 + 3*q1]*dg[k1 + 3*q4] + w*g[q1]*g[q3]*g[q5]*dg[k2 + 3*q2]*dg[k1 + 3*q4] + w*g[q1]*g[q2]*g[q5]*dg[k2 + 3*q3]*dg[k1 + 3*q4] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q5]*dg[k2 + 3*q4] + w*g[q2]*g[q3]*g[q5]*dg[k1 + 3*q1]*dg[k2 + 3*q4] + w*g[q1]*g[q3]*g[q5]*dg[k1 + 3*q2]*dg[k2 + 3*q4] + w*g[q1]*g[q2]*g[q5]*dg[k1 + 3*q3]*dg[k2 + 3*q4] + w*g[q2]*g[q3]*g[q4]*dg[k2 + 3*q1]*dg[k1 + 3*q5] + w*g[q1]*g[q3]*g[q4]*dg[k2 + 3*q2]*dg[k1 + 3*q5] + w*g[q1]*g[q2]*g[q4]*dg[k2 + 3*q3]*dg[k1 + 3*q5] + w*g[q1]*g[q2]*g[q3]*dg[k2 + 3*q4]*dg[k1 + 3*q5] + dw[k2]*(g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + 3*q2] + g[q2]*(g[q1]*g[q4]*g[q5]*dg[k1 + 3*q3] + g[q3]*(g[q4]*g[q5]*dg[k1 + 3*q1] + g[q1]*g[q5]*dg[k1 + 3*q4] + g[q1]*g[q4]*dg[k1 + 3*q5]))) + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*dg[k2 + 3*q5] + w*g[q2]*g[q3]*g[q4]*dg[k1 + 3*q1]*dg[k2 + 3*q5] + w*g[q1]*g[q3]*g[q4]*dg[k1 + 3*q2]*dg[k2 + 3*q5] + w*g[q1]*g[q2]*g[q4]*dg[k1 + 3*q3]*dg[k2 + 3*q5] + w*g[q1]*g[q2]*g[q3]*dg[k1 + 3*q4]*dg[k2 + 3*q5]);
                  for (auto q6 = q5; q6 < dim; ++q6) {
                    ddm6[2187*q1 + 729*q2 + 243*q3 + 81*q4 + 27*q5 + 9*q6 + 3*k1] += v*(ddw[3*k1 + k2]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6] + dw[k1]*g[q2]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + 3*q1] + w*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + 3*q1]*dg[k1 + 3*q2] + dw[k1]*g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k2 + 3*q2] + w*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + 3*q1]*dg[k2 + 3*q2] + w*g[q2]*g[q4]*g[q5]*g[q6]*dg[k2 + 3*q1]*dg[k1 + 3*q3] + w*g[q1]*g[q4]*g[q5]*g[q6]*dg[k2 + 3*q2]*dg[k1 + 3*q3] + dw[k1]*g[q1]*g[q2]*g[q4]*g[q5]*g[q6]*dg[k2 + 3*q3] + w*g[q2]*g[q4]*g[q5]*g[q6]*dg[k1 + 3*q1]*dg[k2 + 3*q3] + w*g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + 3*q2]*dg[k2 + 3*q3] + w*g[q2]*g[q3]*g[q5]*g[q6]*dg[k2 + 3*q1]*dg[k1 + 3*q4] + w*g[q1]*g[q3]*g[q5]*g[q6]*dg[k2 + 3*q2]*dg[k1 + 3*q4] + w*g[q1]*g[q2]*g[q5]*g[q6]*dg[k2 + 3*q3]*dg[k1 + 3*q4] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q5]*g[q6]*dg[k2 + 3*q4] + w*g[q2]*g[q3]*g[q5]*g[q6]*dg[k1 + 3*q1]*dg[k2 + 3*q4] + w*g[q1]*g[q3]*g[q5]*g[q6]*dg[k1 + 3*q2]*dg[k2 + 3*q4] + w*g[q1]*g[q2]*g[q5]*g[q6]*dg[k1 + 3*q3]*dg[k2 + 3*q4] + w*g[q2]*g[q3]*g[q4]*g[q6]*dg[k2 + 3*q1]*dg[k1 + 3*q5] + w*g[q1]*g[q3]*g[q4]*g[q6]*dg[k2 + 3*q2]*dg[k1 + 3*q5] + w*g[q1]*g[q2]*g[q4]*g[q6]*dg[k2 + 3*q3]*dg[k1 + 3*q5] + w*g[q1]*g[q2]*g[q3]*g[q6]*dg[k2 + 3*q4]*dg[k1 + 3*q5] + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q6]*dg[k2 + 3*q5] + w*g[q2]*g[q3]*g[q4]*g[q6]*dg[k1 + 3*q1]*dg[k2 + 3*q5] + w*g[q1]*g[q3]*g[q4]*g[q6]*dg[k1 + 3*q2]*dg[k2 + 3*q5] + w*g[q1]*g[q2]*g[q4]*g[q6]*dg[k1 + 3*q3]*dg[k2 + 3*q5] + w*g[q1]*g[q2]*g[q3]*g[q6]*dg[k1 + 3*q4]*dg[k2 + 3*q5] + w*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + 3*q1]*dg[k1 + 3*q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k2 + 3*q2]*dg[k1 + 3*q6] + w*g[q1]*g[q2]*g[q4]*g[q5]*dg[k2 + 3*q3]*dg[k1 + 3*q6] + w*g[q1]*g[q2]*g[q3]*g[q5]*dg[k2 + 3*q4]*dg[k1 + 3*q6] + w*g[q1]*g[q2]*g[q3]*g[q4]*dg[k2 + 3*q5]*dg[k1 + 3*q6] + dw[k2]*(g[q1]*g[q3]*g[q4]*g[q5]*g[q6]*dg[k1 + 3*q2] + g[q2]*(g[q1]*g[q4]*g[q5]*g[q6]*dg[k1 + 3*q3] + g[q3]*(g[q1]*g[q5]*g[q6]*dg[k1 + 3*q4] + g[q4]*(g[q5]*g[q6]*dg[k1 + 3*q1] + g[q1]*g[q6]*dg[k1 + 3*q5] + g[q1]*g[q5]*dg[k1 + 3*q6])))) + dw[k1]*g[q1]*g[q2]*g[q3]*g[q4]*g[q5]*dg[k2 + 3*q6] + w*g[q2]*g[q3]*g[q4]*g[q5]*dg[k1 + 3*q1]*dg[k2 + 3*q6] + w*g[q1]*g[q3]*g[q4]*g[q5]*dg[k1 + 3*q2]*dg[k2 + 3*q6] + w*g[q1]*g[q2]*g[q4]*g[q5]*dg[k1 + 3*q3]*dg[k2 + 3*q6] + w*g[q1]*g[q2]*g[q3]*g[q5]*dg[k1 + 3*q4]*dg[k2 + 3*q6] + w*g[q1]*g[q2]*g[q3]*g[q4]*dg[k1 + 3*q5]*dg[k2 + 3*q6]);
                  }
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
// Fill in the symmetries for the tensors
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
fillSymmetries(typename Dimension::ThirdRankTensor& ddb) {
  const auto dim = Dimension::nDim;

  // Fill in derivatives that were not assigned
  for (auto q1 = 0; q1 < dim; ++q1) {
    for (auto k1 = 1; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < k1; ++k2) {
        ddb[dim*dim*q1 + dim*k1 + k2] = ddb[dim*dim*q1 + dim*k2 + k1];
      }
    }
  }
}

template<typename Dimension>
inline
void
fillSymmetries(typename Dimension::Tensor& c,
               typename Dimension::ThirdRankTensor& dc,
               typename Dimension::FourthRankTensor& ddc) {
  const auto dim = Dimension::nDim;

  // For q1, q2 that were assigned, fill in derivatives for k1, k2 that were not assigned
  for (auto q1 = 0; q1 < dim; ++q1) {
    for (auto q2 = q1; q2 < dim; ++q2) {
      for (auto k1 = 1; k1 < dim; ++k1) {
        for (auto k2 = 0; k2 < k1; ++k2) {
          ddc[dim*dim*dim*q1 + dim*dim*q2 + dim*k1 + k2] = ddc[dim*dim*dim*q1 + dim*dim*q2 + dim*k2 + k1];
        }
      }
    }
  }
  
  // For q1, q2 that were not assigned, fill in all k1, k2
  for (auto q1 = 1; q1 < dim; ++q1) {
    for (auto q2 = 0; q2 < q1; ++q2) {
      c[dim*q1 + q2] = c[dim*q2 + q1];
      for (auto k1 = 0; k1 < dim; ++k1) {
        dc[dim*dim*q1 + dim*q2 + k1] = dc[dim*dim*q2 + dim*q1 + k1];
        for (auto k2 = 0; k2 < dim; ++k2) {
          ddc[dim*dim*dim*q1 + dim*dim*q2 + dim*k1 + k2] = ddc[dim*dim*dim*q2 + dim*dim*q1 + dim*k1 + k2];
        }
      }
    }
  }
}

template<typename Dimension>
inline
void
fillSymmetries(typename Dimension::ThirdRankTensor& d,
               typename Dimension::FourthRankTensor& dd,
               typename Dimension::FifthRankTensor& ddd) {
  const auto dim = Dimension::nDim;

  // Just sort the indices to get the one that has been filled in
  // This is easy and inefficient
  auto sort2 = [](int& p1, int& p2) {
    if (p1 > p2) {
      std::swap(p1, p2);
    }
    return;
  };
  auto sort3 = [](int& p1, int& p2, int& p3) {
    if (p1 > p2) {
      std::swap(p1, p2);
    }
    if (p1 > p3) {
      std::swap(p1, p3);
    }
    if (p2 > p3) {
      std::swap(p2, p3);
    }
    return;
  };
  
  // Get indices sorted such that q1<=q2<=q3 and k1<=k2
  for (auto q1 = 0; q1 < dim; ++q1) {
    for (auto q2 = 0; q2 < dim; ++q2) {
      for (auto q3 = 0; q3 < dim; ++q3) {
        auto p1 = q1;
        auto p2 = q2;
        auto p3 = q3;
        sort3(p1, p2, p3);
        d[dim*dim*q1 + dim*q2 + q3] = d[dim*dim*p1 + dim*p2 + p3];
        for (auto k1 = 0; k1 < dim; ++k1) {
          dd[dim*dim*dim*q1 + dim*dim*q2 + dim*q3 + k1] = dd[dim*dim*dim*p1 + dim*dim*p2 + dim*p3 + k1];
          for (auto k2 = 0; k2 < dim; ++k2) {
            auto j1 = k1;
            auto j2 = k2;
            sort2(j1, j2);
            ddd[dim*dim*dim*q1 + dim*dim*dim*q2 + dim*dim*q3 + dim*k1 + k2] = ddd[dim*dim*dim*p1 + dim*dim*dim*p2 + dim*dim*p3 + dim*j1 + j2];;
          }
        }
      }
    }
  }
}
  
//------------------------------------------------------------------------------
// Compute the corrections for a single point
// The derivatives depend on the original functions, meaning that a needs to be
// computed before da, da before dda. Symmetry is exploited in solve such that
// it can be filled in at the end, meaning the ddc calculation only uses those
// elements of dc and c that have been computed, before symmetries are filled in
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
computeCorrections(const RKMomentValues<Dimension>& rkMoments,
                   typename Dimension::Scalar& a,
                   typename Dimension::Vector& da,
                   typename Dimension::Tensor& dda);
template<typename Dimension>
inline
void
computeCorrections(const RKMomentValues<Dimension>& rkMoments,
                   typename Dimension::Scalar& a,
                   typename Dimension::Vector& b,
                   typename Dimension::Vector& da,
                   typename Dimension::Tensor& db,
                   typename Dimension::Tensor& dda,
                   typename Dimension::ThirdRankTensor& ddb);
template<typename Dimension>
inline
void
computeCorrections(const RKMomentValues<Dimension>& rkMoments,
                   typename Dimension::Scalar& a,
                   typename Dimension::Vector& b,
                   typename Dimension::Tensor& c,
                   typename Dimension::Vector& da,
                   typename Dimension::Tensor& db,
                   typename Dimension::ThirdRankTensor& dc,
                   typename Dimension::Tensor& dda,
                   typename Dimension::ThirdRankTensor& ddb,
                   typename Dimension::FourthRankTensor& ddc);
template<typename Dimension>
inline
void
computeCorrections(const RKMomentValues<Dimension>& rkMoments,
                   typename Dimension::Scalar& a,
                   typename Dimension::Vector& b,
                   typename Dimension::Tensor& c,
                   typename Dimension::ThirdRankTensor& d,
                   typename Dimension::Vector& da,
                   typename Dimension::Tensor& db,
                   typename Dimension::ThirdRankTensor& dc,
                   typename Dimension::FourthRankTensor& dd,
                   typename Dimension::Tensor& dda,
                   typename Dimension::ThirdRankTensor& ddb,
                   typename Dimension::FourthRankTensor& ddc,
                   typename Dimension::FifthRankTensor& ddd);
// d=1, o=0
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<1>>& rkMoments,
                   Dim<1>::Scalar& a,
                   Dim<1>::Vector& da,
                   Dim<1>::Tensor& dda) {
  const auto dim = Dim<1>::nDim;
  const auto& m0 = rkMoments.m0;
  const auto& dm0 = rkMoments.dm0;
  const auto& ddm0 = rkMoments.ddm0;
  const auto k1 = 0;
  const auto k2 = 0;
  
  // Value
  a = 1./m0;
  // Gradients
  da[0] = -a*dm0[k1]/m0;
  // Hessians
  dda[0] = -(a*ddm0[k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2])/m0;
}

// d=2, o=0
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<2>>& rkMoments,
                   Dim<2>::Scalar& a,
                   Dim<2>::Vector& da,
                   Dim<2>::Tensor& dda) {
  
  const auto dim = Dim<2>::nDim;
  const auto& m0 = rkMoments.m0;
  const auto& dm0 = rkMoments.dm0;
  const auto& ddm0 = rkMoments.ddm0;
  
  // Value
  a = 1./m0;
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    da[k1] = -a*dm0[k1]/m0;
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      dda[2*k1 + k2] = -(a*ddm0[2*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2])/m0;
    }
  }
}
// d=3, o=0
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<3>>& rkMoments,
                   Dim<3>::Scalar& a,
                   Dim<3>::Vector& da,
                   Dim<3>::Tensor& dda) {
  
  const auto dim = Dim<3>::nDim;
  const auto& m0 = rkMoments.m0;
  const auto& dm0 = rkMoments.dm0;
  const auto& ddm0 = rkMoments.ddm0;
  
  // Value
  a = 1./m0;
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    da[k1] = -a*dm0[k1]/m0;
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      dda[3*k1 + k2] = -(a*ddm0[3*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2])/m0;
    }
  }
}
// d=1, o=1
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<1>>& rkMoments,
                   Dim<1>::Scalar& a,
                   Dim<1>::Vector& b,
                   Dim<1>::Vector& da,
                   Dim<1>::Tensor& db,
                   Dim<1>::Tensor& dda,
                   Dim<1>::ThirdRankTensor& ddb) {
  const auto dim = Dim<1>::nDim;
  const auto size = dim + 1;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;
  const auto k1 = 0;
  const auto k2 = 0;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix <<
    m0,    m1[0],
    m1[0], m2[0];
  auto solver = matrix.colPivHouseholderQr();

  // Solve for values
  rhs(0) = 1.;
  rhs(1) = 0.;
  lhs = solver.solve(rhs);
  a = lhs(0);
  b[0] = lhs(1);
  // Solve for derivatives
  rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1]);
  rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1]);
  lhs = solver.solve(rhs);
  da[k1] = lhs(0);
  db[k1] = lhs(1);
  // Solve for second derivatives
  rhs(0) = -(a*ddm0[k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2]);
  rhs(1) = -(a*ddm1[k1 + k2] + b[0]*ddm2[k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2]);
  lhs = solver.solve(rhs);
  dda[k1 + k2] = lhs(0);
  ddb[k1 + k2] = lhs(1);
}
// d=2, o=1
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<2>>& rkMoments,
                   Dim<2>::Scalar& a,
                   Dim<2>::Vector& b,
                   Dim<2>::Vector& da,
                   Dim<2>::Tensor& db,
                   Dim<2>::Tensor& dda,
                   Dim<2>::ThirdRankTensor& ddb) {
  const auto dim = Dim<2>::nDim;
  const auto size = dim + 1;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix <<
    m0, m1[0], m1[1], 
    m1[0], m2[0], m2[1], 
    m1[1], m2[1], m2[3];
  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = 1;
  rhs(1) = 0;
  rhs(2) = 0;
  
  lhs = solver.solve(rhs);
  
  a = lhs(0);
  b[0] = lhs(1);
  b[1] = lhs(2);
  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[2 + k1]);
    rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + b[1]*dm2[4 + k1]);
    rhs(2) = -(a*dm1[2 + k1] + b[0]*dm2[2 + k1] + b[1]*dm2[6 + k1]);
    
    lhs = solver.solve(rhs);
    
    da[k1] = lhs(0);
    db[k1] = lhs(1);
    db[2 + k1] = lhs(2);
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[2*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[2*k1 + k2] + b[1]*ddm1[4 + 2*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[2 + k2]*dm1[2 + k1] + db[2 + k1]*dm1[2 + k2]);
      rhs(1) = -(a*ddm1[2*k1 + k2] + b[0]*ddm2[2*k1 + k2] + b[1]*ddm2[8 + 2*k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + db[2 + k2]*dm2[4 + k1] + db[2 + k1]*dm2[4 + k2]);
      rhs(2) = -(a*ddm1[4 + 2*k1 + k2] + b[0]*ddm2[4 + 2*k1 + k2] + b[1]*ddm2[12 + 2*k1 + k2] + da[k2]*dm1[2 + k1] + da[k1]*dm1[2 + k2] + db[k2]*dm2[2 + k1] + db[k1]*dm2[2 + k2] + db[2 + k2]*dm2[6 + k1] + db[2 + k1]*dm2[6 + k2]);
      
      lhs = solver.solve(rhs);
      
      dda[2*k1 + k2] = lhs(0);
      ddb[2*k1 + k2] = lhs(1);
      ddb[4 + 2*k1 + k2] = lhs(2);
    }
  }

  // Fill symmetries
  fillSymmetries<Dim<2>>(ddb);
}
// d=3, o=1
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<3>>& rkMoments,
                   Dim<3>::Scalar& a,
                   Dim<3>::Vector& b,
                   Dim<3>::Vector& da,
                   Dim<3>::Tensor& db,
                   Dim<3>::Tensor& dda,
                   Dim<3>::ThirdRankTensor& ddb) {
  const auto dim = Dim<3>::nDim;
  const auto size = dim + 1;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix <<
    m0, m1[0], m1[1], m1[2], 
    m1[0], m2[0], m2[1], m2[2], 
    m1[1], m2[1], m2[4], m2[5], 
    m1[2], m2[2], m2[5], m2[8];

  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = (1);
  rhs(1) = (0);
  rhs(2) = (0);
  rhs(3) = (0);
  
  lhs = solver.solve(rhs);

  a = lhs(0);
  b[0] = lhs(1);
  b[1] = lhs(2);
  b[2] = lhs(3);
  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[3 + k1] + b[2]*dm1[6 + k1]);
    rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + b[1]*dm2[9 + k1] + b[2]*dm2[18 + k1]);
    rhs(2) = -(a*dm1[3 + k1] + b[0]*dm2[3 + k1] + b[1]*dm2[12 + k1] + b[2]*dm2[21 + k1]);
    rhs(3) = -(a*dm1[6 + k1] + b[0]*dm2[6 + k1] + b[1]*dm2[15 + k1] + b[2]*dm2[24 + k1]);

    lhs = solver.solve(rhs);

    da[k1] = lhs(0);
    db[k1] = lhs(1);
    db[3 + k1] = lhs(2);
    db[6 + k1] = lhs(3);
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[3*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[3*k1 + k2] + b[1]*ddm1[9 + 3*k1 + k2] + b[2]*ddm1[18 + 3*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[3 + k2]*dm1[3 + k1] + db[3 + k1]*dm1[3 + k2] + db[6 + k2]*dm1[6 + k1] + db[6 + k1]*dm1[6 + k2]);
      rhs(1) = -(a*ddm1[3*k1 + k2] + b[0]*ddm2[3*k1 + k2] + b[1]*ddm2[27 + 3*k1 + k2] + b[2]*ddm2[54 + 3*k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + db[3 + k2]*dm2[9 + k1] + db[3 + k1]*dm2[9 + k2] + db[6 + k2]*dm2[18 + k1] + db[6 + k1]*dm2[18 + k2]);
      rhs(2) = -(a*ddm1[9 + 3*k1 + k2] + b[0]*ddm2[9 + 3*k1 + k2] + b[1]*ddm2[36 + 3*k1 + k2] + b[2]*ddm2[63 + 3*k1 + k2] + da[k2]*dm1[3 + k1] + da[k1]*dm1[3 + k2] + db[k2]*dm2[3 + k1] + db[k1]*dm2[3 + k2] + db[3 + k2]*dm2[12 + k1] + db[3 + k1]*dm2[12 + k2] + db[6 + k2]*dm2[21 + k1] + db[6 + k1]*dm2[21 + k2]);
      rhs(3) = -(a*ddm1[18 + 3*k1 + k2] + b[0]*ddm2[18 + 3*k1 + k2] + b[1]*ddm2[45 + 3*k1 + k2] + b[2]*ddm2[72 + 3*k1 + k2] + da[k2]*dm1[6 + k1] + da[k1]*dm1[6 + k2] + db[k2]*dm2[6 + k1] + db[k1]*dm2[6 + k2] + db[3 + k2]*dm2[15 + k1] + db[3 + k1]*dm2[15 + k2] + db[6 + k2]*dm2[24 + k1] + db[6 + k1]*dm2[24 + k2]);
      
      lhs = solver.solve(rhs);

      dda[3*k1 + k2] = lhs(0);
      ddb[3*k1 + k2] = lhs(1);
      ddb[9 + 3*k1 + k2] = lhs(2);
      ddb[18 + 3*k1 + k2] = lhs(3);      
    }
  }

  // Fill symmetries
  fillSymmetries<Dim<3>>(ddb);
}
// d=1, o=2
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<1>>& rkMoments,
                   Dim<1>::Scalar& a,
                   Dim<1>::Vector& b,
                   Dim<1>::Tensor& c,
                   Dim<1>::Vector& da,
                   Dim<1>::Tensor& db,
                   Dim<1>::ThirdRankTensor& dc,
                   Dim<1>::Tensor& dda,
                   Dim<1>::ThirdRankTensor& ddb,
                   Dim<1>::FourthRankTensor& ddc) {
  const auto dim = Dim<1>::nDim;
  const auto size = 3;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& m3 = rkMoments.m3;
  const auto& m4 = rkMoments.m4;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& dm3 = rkMoments.dm3;
  const auto& dm4 = rkMoments.dm4;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;
  const auto& ddm3 = rkMoments.ddm3;
  const auto& ddm4 = rkMoments.ddm4;
  const auto k1 = 0;
  const auto k2 = 0;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix << 
    m0,    m1[0], m2[0],
    m1[0], m2[0], m3[0],
    m2[0], m3[0], m4[0];
  
  auto solver = matrix.colPivHouseholderQr();

  // Solve for values
  rhs(0) = 1.;
  rhs(1) = 0.;
  rhs(2) = 0.;
  
  lhs = solver.solve(rhs);
  
  a = lhs(0);
  b[0] = lhs(1);
  c[0] = lhs(2);
  // Solve for derivatives
  rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + c[0]*dm2[k1]);
  rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + c[0]*dm3[k1]);
  rhs(2) = -(a*dm2[k1] + b[0]*dm3[k1] + c[0]*dm4[k1]);
  
  lhs = solver.solve(rhs);
  
  da[0] = lhs(0);
  db[0] = lhs(1);
  dc[0] = lhs(2);
  // Solve for second derivatives
  rhs(0) = -(a*ddm0[k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[k1 + k2] + c[0]*ddm2[k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + dc[k2]*dm2[k1] + dc[k1]*dm2[k2]);
  rhs(1) = -(a*ddm1[k1 + k2] + b[0]*ddm2[k1 + k2] + c[0]*ddm3[k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + dc[k2]*dm3[k1] + dc[k1]*dm3[k2]);
  rhs(2) = -(a*ddm2[k1 + k2] + b[0]*ddm3[k1 + k2] + c[0]*ddm4[k1 + k2] + da[k2]*dm2[k1] + da[k1]*dm2[k2] + db[k2]*dm3[k1] + db[k1]*dm3[k2] + dc[k2]*dm4[k1] + dc[k1]*dm4[k2]);

  lhs = solver.solve(rhs);

  dda[k1 + k2] = lhs(0);
  ddb[k1 + k2] = lhs(1);
  ddc[k1 + k2] = lhs(2);
}
// d=2, o=2
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<2>>& rkMoments,
                   Dim<2>::Scalar& a,
                   Dim<2>::Vector& b,
                   Dim<2>::Tensor& c,
                   Dim<2>::Vector& da,
                   Dim<2>::Tensor& db,
                   Dim<2>::ThirdRankTensor& dc,
                   Dim<2>::Tensor& dda,
                   Dim<2>::ThirdRankTensor& ddb,
                   Dim<2>::FourthRankTensor& ddc) {
  const auto dim = Dim<2>::nDim;
  const auto size = 6;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& m3 = rkMoments.m3;
  const auto& m4 = rkMoments.m4;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& dm3 = rkMoments.dm3;
  const auto& dm4 = rkMoments.dm4;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;
  const auto& ddm3 = rkMoments.ddm3;
  const auto& ddm4 = rkMoments.ddm4;

  // Get matrix to invert
  VectorType rhs;
  MatrixType matrix;
  VectorType lhs;
  matrix <<
    m0, m1[0], m1[1], m2[0], 2*m2[1], m2[3], 
    m1[0], m2[0], m2[1], m3[0], 2*m3[1], m3[3], 
    m1[1], m2[1], m2[3], m3[1], 2*m3[3], m3[7], 
    m2[0], m3[0], m3[1], m4[0], 2*m4[1], m4[3], 
    m2[1], m3[1], m3[3], m4[1], 2*m4[3], m4[7], 
    m2[3], m3[3], m3[7], m4[3], 2*m4[7], m4[15];
  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = 1;
  rhs(1) = 0;
  rhs(2) = 0;
  rhs(3) = 0;
  rhs(4) = 0;
  rhs(5) = 0;
 
  lhs = solver.solve(rhs);
  a = lhs(0);
  b[0] = lhs(1);
  b[1] = lhs(2);
  c[0] = lhs(3);
  c[1] = lhs(4);
  c[3] = lhs(5);
  
  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[2 + k1] + c[0]*dm2[k1] + 2*c[1]*dm2[2 + k1] + c[3]*dm2[6 + k1]);
    rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + b[1]*dm2[2 + k1] + c[0]*dm3[k1] + 2*c[1]*dm3[2 + k1] + c[3]*dm3[6 + k1]);
    rhs(2) = -(a*dm1[2 + k1] + b[0]*dm2[2 + k1] + b[1]*dm2[6 + k1] + c[0]*dm3[2 + k1] + 2*c[1]*dm3[6 + k1] + c[3]*dm3[14 + k1]);
    rhs(3) = -(a*dm2[k1] + b[0]*dm3[k1] + b[1]*dm3[2 + k1] + c[0]*dm4[k1] + 2*c[1]*dm4[2 + k1] + c[3]*dm4[6 + k1]);
    rhs(4) = -(a*dm2[2 + k1] + b[0]*dm3[2 + k1] + b[1]*dm3[6 + k1] + c[0]*dm4[2 + k1] + 2*c[1]*dm4[6 + k1] + c[3]*dm4[14 + k1]);
    rhs(5) = -(a*dm2[6 + k1] + b[0]*dm3[6 + k1] + b[1]*dm3[14 + k1] + c[0]*dm4[6 + k1] + 2*c[1]*dm4[14 + k1] + c[3]*dm4[30 + k1]);
 
    lhs = solver.solve(rhs);
      
    da[k1] = lhs(0);
    db[k1] = lhs(1);
    db[2 + k1] = lhs(2);
    dc[k1] = lhs(3);
    dc[2 + k1] = lhs(4);
    dc[6 + k1] = lhs(5);
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[2*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[2*k1 + k2] + b[1]*ddm1[4 + 2*k1 + k2] + c[0]*ddm2[2*k1 + k2] + 2*c[1]*ddm2[4 + 2*k1 + k2] + c[3]*ddm2[12 + 2*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[2 + k2]*dm1[2 + k1] + db[2 + k1]*dm1[2 + k2] + dc[k2]*dm2[k1] + dc[k1]*dm2[k2] + 2*dc[2 + k2]*dm2[2 + k1] + 2*dc[2 + k1]*dm2[2 + k2] + dc[6 + k2]*dm2[6 + k1] + dc[6 + k1]*dm2[6 + k2]);
      rhs(1) = -(a*ddm1[2*k1 + k2] + b[0]*ddm2[2*k1 + k2] + b[1]*ddm2[4 + 2*k1 + k2] + c[0]*ddm3[2*k1 + k2] + 2*c[1]*ddm3[4 + 2*k1 + k2] + c[3]*ddm3[12 + 2*k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + db[2 + k2]*dm2[2 + k1] + db[2 + k1]*dm2[2 + k2] + dc[k2]*dm3[k1] + dc[k1]*dm3[k2] + 2*dc[2 + k2]*dm3[2 + k1] + 2*dc[2 + k1]*dm3[2 + k2] + dc[6 + k2]*dm3[6 + k1] + dc[6 + k1]*dm3[6 + k2]);
      rhs(2) = -(a*ddm1[4 + 2*k1 + k2] + b[0]*ddm2[4 + 2*k1 + k2] + b[1]*ddm2[12 + 2*k1 + k2] + c[0]*ddm3[4 + 2*k1 + k2] + 2*c[1]*ddm3[12 + 2*k1 + k2] + c[3]*ddm3[28 + 2*k1 + k2] + da[k2]*dm1[2 + k1] + da[k1]*dm1[2 + k2] + db[k2]*dm2[2 + k1] + db[k1]*dm2[2 + k2] + db[2 + k2]*dm2[6 + k1] + db[2 + k1]*dm2[6 + k2] + dc[k2]*dm3[2 + k1] + dc[k1]*dm3[2 + k2] + 2*dc[2 + k2]*dm3[6 + k1] + 2*dc[2 + k1]*dm3[6 + k2] + dc[6 + k2]*dm3[14 + k1] + dc[6 + k1]*dm3[14 + k2]);
      rhs(3) = -(a*ddm2[2*k1 + k2] + b[0]*ddm3[2*k1 + k2] + b[1]*ddm3[4 + 2*k1 + k2] + c[0]*ddm4[2*k1 + k2] + 2*c[1]*ddm4[4 + 2*k1 + k2] + c[3]*ddm4[12 + 2*k1 + k2] + da[k2]*dm2[k1] + da[k1]*dm2[k2] + db[k2]*dm3[k1] + db[k1]*dm3[k2] + db[2 + k2]*dm3[2 + k1] + db[2 + k1]*dm3[2 + k2] + dc[k2]*dm4[k1] + dc[k1]*dm4[k2] + 2*dc[2 + k2]*dm4[2 + k1] + 2*dc[2 + k1]*dm4[2 + k2] + dc[6 + k2]*dm4[6 + k1] + dc[6 + k1]*dm4[6 + k2]);
      rhs(4) = -(a*ddm2[4 + 2*k1 + k2] + b[0]*ddm3[4 + 2*k1 + k2] + b[1]*ddm3[12 + 2*k1 + k2] + c[0]*ddm4[4 + 2*k1 + k2] + 2*c[1]*ddm4[12 + 2*k1 + k2] + c[3]*ddm4[28 + 2*k1 + k2] + da[k2]*dm2[2 + k1] + da[k1]*dm2[2 + k2] + db[k2]*dm3[2 + k1] + db[k1]*dm3[2 + k2] + db[2 + k2]*dm3[6 + k1] + db[2 + k1]*dm3[6 + k2] + dc[k2]*dm4[2 + k1] + dc[k1]*dm4[2 + k2] + 2*dc[2 + k2]*dm4[6 + k1] + 2*dc[2 + k1]*dm4[6 + k2] + dc[6 + k2]*dm4[14 + k1] + dc[6 + k1]*dm4[14 + k2]);
      rhs(5) = -(a*ddm2[12 + 2*k1 + k2] + b[0]*ddm3[12 + 2*k1 + k2] + b[1]*ddm3[28 + 2*k1 + k2] + c[0]*ddm4[12 + 2*k1 + k2] + 2*c[1]*ddm4[28 + 2*k1 + k2] + c[3]*ddm4[60 + 2*k1 + k2] + da[k2]*dm2[6 + k1] + da[k1]*dm2[6 + k2] + db[k2]*dm3[6 + k1] + db[k1]*dm3[6 + k2] + db[2 + k2]*dm3[14 + k1] + db[2 + k1]*dm3[14 + k2] + dc[k2]*dm4[6 + k1] + dc[k1]*dm4[6 + k2] + 2*dc[2 + k2]*dm4[14 + k1] + 2*dc[2 + k1]*dm4[14 + k2] + dc[6 + k2]*dm4[30 + k1] + dc[6 + k1]*dm4[30 + k2]);

      lhs = solver.solve(rhs);

      dda[2*k1 + k2] = lhs(0);
      ddb[2*k1 + k2] = lhs(1);
      ddb[4 + 2*k1 + k2] = lhs(2);
      ddc[2*k1 + k2] = lhs(3);
      ddc[4 + 2*k1 + k2] = lhs(4);
      ddc[12 + 2*k1 + k2] = lhs(5);
    }
  }

  // Fill symmetries
  fillSymmetries<Dim<2>>(ddb);
  fillSymmetries<Dim<2>>(c, dc, ddc);
}
// d=3, o=2
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<3>>& rkMoments,
                   Dim<3>::Scalar& a,
                   Dim<3>::Vector& b,
                   Dim<3>::Tensor& c,
                   Dim<3>::Vector& da,
                   Dim<3>::Tensor& db,
                   Dim<3>::ThirdRankTensor& dc,
                   Dim<3>::Tensor& dda,
                   Dim<3>::ThirdRankTensor& ddb,
                   Dim<3>::FourthRankTensor& ddc) {
  const auto dim = Dim<3>::nDim;
  const auto size = 10;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& m3 = rkMoments.m3;
  const auto& m4 = rkMoments.m4;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& dm3 = rkMoments.dm3;
  const auto& dm4 = rkMoments.dm4;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;
  const auto& ddm3 = rkMoments.ddm3;
  const auto& ddm4 = rkMoments.ddm4;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix <<
    m0, m1[0], m1[1], m1[2], m2[0], 2*m2[1], 2*m2[2], m2[4], 2*m2[5], m2[8], 
    m1[0], m2[0], m2[1], m2[2], m3[0], 2*m3[1], 2*m3[2], m3[4], 2*m3[5], m3[8], 
    m1[1], m2[1], m2[4], m2[5], m3[1], 2*m3[4], 2*m3[5], m3[13], 2*m3[14], m3[17], 
    m1[2], m2[2], m2[5], m2[8], m3[2], 2*m3[5], 2*m3[8], m3[14], 2*m3[17], m3[26], 
    m2[0], m3[0], m3[1], m3[2], m4[0], 2*m4[1], 2*m4[2], m4[4], 2*m4[5], m4[8], 
    m2[1], m3[1], m3[4], m3[5], m4[1], 2*m4[4], 2*m4[5], m4[13], 2*m4[14], m4[17], 
    m2[2], m3[2], m3[5], m3[8], m4[2], 2*m4[5], 2*m4[8], m4[14], 2*m4[17], m4[26], 
    m2[4], m3[4], m3[13], m3[14], m4[4], 2*m4[13], 2*m4[14], m4[40], 2*m4[41], m4[44], 
    m2[5], m3[5], m3[14], m3[17], m4[5], 2*m4[14], 2*m4[17], m4[41], 2*m4[44], m4[53], 
    m2[8], m3[8], m3[17], m3[26], m4[8], 2*m4[17], 2*m4[26], m4[44], 2*m4[53], m4[80];
  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = 1;
  rhs(1) = 0;
  rhs(2) = 0;
  rhs(3) = 0;
  rhs(4) = 0;
  rhs(5) = 0;
  rhs(6) = 0;
  rhs(7) = 0;
  rhs(8) = 0;
  rhs(9) = 0;
 
  lhs = solver.solve(rhs);

  a = lhs(0);
  b[0] = lhs(1);
  b[1] = lhs(2);
  b[2] = lhs(3);
  c[0] = lhs(4);
  c[1] = lhs(5);
  c[2] = lhs(6);
  c[4] = lhs(7);
  c[5] = lhs(8);
  c[8] = lhs(9);  
  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[3 + k1] + b[2]*dm1[6 + k1] + c[0]*dm2[k1] + 2*c[1]*dm2[3 + k1] + 2*c[2]*dm2[6 + k1] + c[4]*dm2[12 + k1] + 2*c[5]*dm2[15 + k1] + c[8]*dm2[24 + k1]);
    rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + b[1]*dm2[3 + k1] + b[2]*dm2[6 + k1] + c[0]*dm3[k1] + 2*c[1]*dm3[3 + k1] + 2*c[2]*dm3[6 + k1] + c[4]*dm3[12 + k1] + 2*c[5]*dm3[15 + k1] + c[8]*dm3[24 + k1]);
    rhs(2) = -(a*dm1[3 + k1] + b[0]*dm2[3 + k1] + b[1]*dm2[12 + k1] + b[2]*dm2[15 + k1] + c[0]*dm3[3 + k1] + 2*c[1]*dm3[12 + k1] + 2*c[2]*dm3[15 + k1] + c[4]*dm3[39 + k1] + 2*c[5]*dm3[42 + k1] + c[8]*dm3[51 + k1]);
    rhs(3) = -(a*dm1[6 + k1] + b[0]*dm2[6 + k1] + b[1]*dm2[15 + k1] + b[2]*dm2[24 + k1] + c[0]*dm3[6 + k1] + 2*c[1]*dm3[15 + k1] + 2*c[2]*dm3[24 + k1] + c[4]*dm3[42 + k1] + 2*c[5]*dm3[51 + k1] + c[8]*dm3[78 + k1]);
    rhs(4) = -(a*dm2[k1] + b[0]*dm3[k1] + b[1]*dm3[3 + k1] + b[2]*dm3[6 + k1] + c[0]*dm4[k1] + 2*c[1]*dm4[3 + k1] + 2*c[2]*dm4[6 + k1] + c[4]*dm4[12 + k1] + 2*c[5]*dm4[15 + k1] + c[8]*dm4[24 + k1]);
    rhs(5) = -(a*dm2[3 + k1] + b[0]*dm3[3 + k1] + b[1]*dm3[12 + k1] + b[2]*dm3[15 + k1] + c[0]*dm4[3 + k1] + 2*c[1]*dm4[12 + k1] + 2*c[2]*dm4[15 + k1] + c[4]*dm4[39 + k1] + 2*c[5]*dm4[42 + k1] + c[8]*dm4[51 + k1]);
    rhs(6) = -(a*dm2[6 + k1] + b[0]*dm3[6 + k1] + b[1]*dm3[15 + k1] + b[2]*dm3[24 + k1] + c[0]*dm4[6 + k1] + 2*c[1]*dm4[15 + k1] + 2*c[2]*dm4[24 + k1] + c[4]*dm4[42 + k1] + 2*c[5]*dm4[51 + k1] + c[8]*dm4[78 + k1]);
    rhs(7) = -(a*dm2[12 + k1] + b[0]*dm3[12 + k1] + b[1]*dm3[39 + k1] + b[2]*dm3[42 + k1] + c[0]*dm4[12 + k1] + 2*c[1]*dm4[39 + k1] + 2*c[2]*dm4[42 + k1] + c[4]*dm4[120 + k1] + 2*c[5]*dm4[123 + k1] + c[8]*dm4[132 + k1]);
    rhs(8) = -(a*dm2[15 + k1] + b[0]*dm3[15 + k1] + b[1]*dm3[42 + k1] + b[2]*dm3[51 + k1] + c[0]*dm4[15 + k1] + 2*c[1]*dm4[42 + k1] + 2*c[2]*dm4[51 + k1] + c[4]*dm4[123 + k1] + 2*c[5]*dm4[132 + k1] + c[8]*dm4[159 + k1]);
    rhs(9) = -(a*dm2[24 + k1] + b[0]*dm3[24 + k1] + b[1]*dm3[51 + k1] + b[2]*dm3[78 + k1] + c[0]*dm4[24 + k1] + 2*c[1]*dm4[51 + k1] + 2*c[2]*dm4[78 + k1] + c[4]*dm4[132 + k1] + 2*c[5]*dm4[159 + k1] + c[8]*dm4[240 + k1]);
 
    lhs = solver.solve(rhs);

    da[k1] = lhs(0);
    db[k1] = lhs(1);
    db[3 + k1] = lhs(2);
    db[6 + k1] = lhs(3);
    dc[k1] = lhs(4);
    dc[3 + k1] = lhs(5);
    dc[6 + k1] = lhs(6);
    dc[12 + k1] = lhs(7);
    dc[15 + k1] = lhs(8);
    dc[24 + k1] = lhs(9);
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[3*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[3*k1 + k2] + b[1]*ddm1[9 + 3*k1 + k2] + b[2]*ddm1[18 + 3*k1 + k2] + c[0]*ddm2[3*k1 + k2] + 2*c[1]*ddm2[9 + 3*k1 + k2] + 2*c[2]*ddm2[18 + 3*k1 + k2] + c[4]*ddm2[36 + 3*k1 + k2] + 2*c[5]*ddm2[45 + 3*k1 + k2] + c[8]*ddm2[72 + 3*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[3 + k2]*dm1[3 + k1] + db[3 + k1]*dm1[3 + k2] + db[6 + k2]*dm1[6 + k1] + db[6 + k1]*dm1[6 + k2] + dc[k2]*dm2[k1] + dc[k1]*dm2[k2] + 2*dc[3 + k2]*dm2[3 + k1] + 2*dc[3 + k1]*dm2[3 + k2] + 2*dc[6 + k2]*dm2[6 + k1] + 2*dc[6 + k1]*dm2[6 + k2] + dc[12 + k2]*dm2[12 + k1] + dc[12 + k1]*dm2[12 + k2] + 2*dc[15 + k2]*dm2[15 + k1] + 2*dc[15 + k1]*dm2[15 + k2] + dc[24 + k2]*dm2[24 + k1] + dc[24 + k1]*dm2[24 + k2]);
      rhs(1) = -(a*ddm1[3*k1 + k2] + b[0]*ddm2[3*k1 + k2] + b[1]*ddm2[9 + 3*k1 + k2] + b[2]*ddm2[18 + 3*k1 + k2] + c[0]*ddm3[3*k1 + k2] + 2*c[1]*ddm3[9 + 3*k1 + k2] + 2*c[2]*ddm3[18 + 3*k1 + k2] + c[4]*ddm3[36 + 3*k1 + k2] + 2*c[5]*ddm3[45 + 3*k1 + k2] + c[8]*ddm3[72 + 3*k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + db[3 + k2]*dm2[3 + k1] + db[3 + k1]*dm2[3 + k2] + db[6 + k2]*dm2[6 + k1] + db[6 + k1]*dm2[6 + k2] + dc[k2]*dm3[k1] + dc[k1]*dm3[k2] + 2*dc[3 + k2]*dm3[3 + k1] + 2*dc[3 + k1]*dm3[3 + k2] + 2*dc[6 + k2]*dm3[6 + k1] + 2*dc[6 + k1]*dm3[6 + k2] + dc[12 + k2]*dm3[12 + k1] + dc[12 + k1]*dm3[12 + k2] + 2*dc[15 + k2]*dm3[15 + k1] + 2*dc[15 + k1]*dm3[15 + k2] + dc[24 + k2]*dm3[24 + k1] + dc[24 + k1]*dm3[24 + k2]);
      rhs(2) = -(a*ddm1[9 + 3*k1 + k2] + b[0]*ddm2[9 + 3*k1 + k2] + b[1]*ddm2[36 + 3*k1 + k2] + b[2]*ddm2[45 + 3*k1 + k2] + c[0]*ddm3[9 + 3*k1 + k2] + 2*c[1]*ddm3[36 + 3*k1 + k2] + 2*c[2]*ddm3[45 + 3*k1 + k2] + c[4]*ddm3[117 + 3*k1 + k2] + 2*c[5]*ddm3[126 + 3*k1 + k2] + c[8]*ddm3[153 + 3*k1 + k2] + da[k2]*dm1[3 + k1] + da[k1]*dm1[3 + k2] + db[k2]*dm2[3 + k1] + db[k1]*dm2[3 + k2] + db[3 + k2]*dm2[12 + k1] + db[3 + k1]*dm2[12 + k2] + db[6 + k2]*dm2[15 + k1] + db[6 + k1]*dm2[15 + k2] + dc[k2]*dm3[3 + k1] + dc[k1]*dm3[3 + k2] + 2*dc[3 + k2]*dm3[12 + k1] + 2*dc[3 + k1]*dm3[12 + k2] + 2*dc[6 + k2]*dm3[15 + k1] + 2*dc[6 + k1]*dm3[15 + k2] + dc[12 + k2]*dm3[39 + k1] + dc[12 + k1]*dm3[39 + k2] + 2*dc[15 + k2]*dm3[42 + k1] + 2*dc[15 + k1]*dm3[42 + k2] + dc[24 + k2]*dm3[51 + k1] + dc[24 + k1]*dm3[51 + k2]);
      rhs(3) = -(a*ddm1[18 + 3*k1 + k2] + b[0]*ddm2[18 + 3*k1 + k2] + b[1]*ddm2[45 + 3*k1 + k2] + b[2]*ddm2[72 + 3*k1 + k2] + c[0]*ddm3[18 + 3*k1 + k2] + 2*c[1]*ddm3[45 + 3*k1 + k2] + 2*c[2]*ddm3[72 + 3*k1 + k2] + c[4]*ddm3[126 + 3*k1 + k2] + 2*c[5]*ddm3[153 + 3*k1 + k2] + c[8]*ddm3[234 + 3*k1 + k2] + da[k2]*dm1[6 + k1] + da[k1]*dm1[6 + k2] + db[k2]*dm2[6 + k1] + db[k1]*dm2[6 + k2] + db[3 + k2]*dm2[15 + k1] + db[3 + k1]*dm2[15 + k2] + db[6 + k2]*dm2[24 + k1] + db[6 + k1]*dm2[24 + k2] + dc[k2]*dm3[6 + k1] + dc[k1]*dm3[6 + k2] + 2*dc[3 + k2]*dm3[15 + k1] + 2*dc[3 + k1]*dm3[15 + k2] + 2*dc[6 + k2]*dm3[24 + k1] + 2*dc[6 + k1]*dm3[24 + k2] + dc[12 + k2]*dm3[42 + k1] + dc[12 + k1]*dm3[42 + k2] + 2*dc[15 + k2]*dm3[51 + k1] + 2*dc[15 + k1]*dm3[51 + k2] + dc[24 + k2]*dm3[78 + k1] + dc[24 + k1]*dm3[78 + k2]);
      rhs(4) = -(a*ddm2[3*k1 + k2] + b[0]*ddm3[3*k1 + k2] + b[1]*ddm3[9 + 3*k1 + k2] + b[2]*ddm3[18 + 3*k1 + k2] + c[0]*ddm4[3*k1 + k2] + 2*c[1]*ddm4[9 + 3*k1 + k2] + 2*c[2]*ddm4[18 + 3*k1 + k2] + c[4]*ddm4[36 + 3*k1 + k2] + 2*c[5]*ddm4[45 + 3*k1 + k2] + c[8]*ddm4[72 + 3*k1 + k2] + da[k2]*dm2[k1] + da[k1]*dm2[k2] + db[k2]*dm3[k1] + db[k1]*dm3[k2] + db[3 + k2]*dm3[3 + k1] + db[3 + k1]*dm3[3 + k2] + db[6 + k2]*dm3[6 + k1] + db[6 + k1]*dm3[6 + k2] + dc[k2]*dm4[k1] + dc[k1]*dm4[k2] + 2*dc[3 + k2]*dm4[3 + k1] + 2*dc[3 + k1]*dm4[3 + k2] + 2*dc[6 + k2]*dm4[6 + k1] + 2*dc[6 + k1]*dm4[6 + k2] + dc[12 + k2]*dm4[12 + k1] + dc[12 + k1]*dm4[12 + k2] + 2*dc[15 + k2]*dm4[15 + k1] + 2*dc[15 + k1]*dm4[15 + k2] + dc[24 + k2]*dm4[24 + k1] + dc[24 + k1]*dm4[24 + k2]);
      rhs(5) = -(a*ddm2[9 + 3*k1 + k2] + b[0]*ddm3[9 + 3*k1 + k2] + b[1]*ddm3[36 + 3*k1 + k2] + b[2]*ddm3[45 + 3*k1 + k2] + c[0]*ddm4[9 + 3*k1 + k2] + 2*c[1]*ddm4[36 + 3*k1 + k2] + 2*c[2]*ddm4[45 + 3*k1 + k2] + c[4]*ddm4[117 + 3*k1 + k2] + 2*c[5]*ddm4[126 + 3*k1 + k2] + c[8]*ddm4[153 + 3*k1 + k2] + da[k2]*dm2[3 + k1] + da[k1]*dm2[3 + k2] + db[k2]*dm3[3 + k1] + db[k1]*dm3[3 + k2] + db[3 + k2]*dm3[12 + k1] + db[3 + k1]*dm3[12 + k2] + db[6 + k2]*dm3[15 + k1] + db[6 + k1]*dm3[15 + k2] + dc[k2]*dm4[3 + k1] + dc[k1]*dm4[3 + k2] + 2*dc[3 + k2]*dm4[12 + k1] + 2*dc[3 + k1]*dm4[12 + k2] + 2*dc[6 + k2]*dm4[15 + k1] + 2*dc[6 + k1]*dm4[15 + k2] + dc[12 + k2]*dm4[39 + k1] + dc[12 + k1]*dm4[39 + k2] + 2*dc[15 + k2]*dm4[42 + k1] + 2*dc[15 + k1]*dm4[42 + k2] + dc[24 + k2]*dm4[51 + k1] + dc[24 + k1]*dm4[51 + k2]);
      rhs(6) = -(a*ddm2[18 + 3*k1 + k2] + b[0]*ddm3[18 + 3*k1 + k2] + b[1]*ddm3[45 + 3*k1 + k2] + b[2]*ddm3[72 + 3*k1 + k2] + c[0]*ddm4[18 + 3*k1 + k2] + 2*c[1]*ddm4[45 + 3*k1 + k2] + 2*c[2]*ddm4[72 + 3*k1 + k2] + c[4]*ddm4[126 + 3*k1 + k2] + 2*c[5]*ddm4[153 + 3*k1 + k2] + c[8]*ddm4[234 + 3*k1 + k2] + da[k2]*dm2[6 + k1] + da[k1]*dm2[6 + k2] + db[k2]*dm3[6 + k1] + db[k1]*dm3[6 + k2] + db[3 + k2]*dm3[15 + k1] + db[3 + k1]*dm3[15 + k2] + db[6 + k2]*dm3[24 + k1] + db[6 + k1]*dm3[24 + k2] + dc[k2]*dm4[6 + k1] + dc[k1]*dm4[6 + k2] + 2*dc[3 + k2]*dm4[15 + k1] + 2*dc[3 + k1]*dm4[15 + k2] + 2*dc[6 + k2]*dm4[24 + k1] + 2*dc[6 + k1]*dm4[24 + k2] + dc[12 + k2]*dm4[42 + k1] + dc[12 + k1]*dm4[42 + k2] + 2*dc[15 + k2]*dm4[51 + k1] + 2*dc[15 + k1]*dm4[51 + k2] + dc[24 + k2]*dm4[78 + k1] + dc[24 + k1]*dm4[78 + k2]);
      rhs(7) = -(a*ddm2[36 + 3*k1 + k2] + b[0]*ddm3[36 + 3*k1 + k2] + b[1]*ddm3[117 + 3*k1 + k2] + b[2]*ddm3[126 + 3*k1 + k2] + c[0]*ddm4[36 + 3*k1 + k2] + 2*c[1]*ddm4[117 + 3*k1 + k2] + 2*c[2]*ddm4[126 + 3*k1 + k2] + c[4]*ddm4[360 + 3*k1 + k2] + 2*c[5]*ddm4[369 + 3*k1 + k2] + c[8]*ddm4[396 + 3*k1 + k2] + da[k2]*dm2[12 + k1] + da[k1]*dm2[12 + k2] + db[k2]*dm3[12 + k1] + db[k1]*dm3[12 + k2] + db[3 + k2]*dm3[39 + k1] + db[3 + k1]*dm3[39 + k2] + db[6 + k2]*dm3[42 + k1] + db[6 + k1]*dm3[42 + k2] + dc[k2]*dm4[12 + k1] + dc[k1]*dm4[12 + k2] + 2*dc[3 + k2]*dm4[39 + k1] + 2*dc[3 + k1]*dm4[39 + k2] + 2*dc[6 + k2]*dm4[42 + k1] + 2*dc[6 + k1]*dm4[42 + k2] + dc[12 + k2]*dm4[120 + k1] + dc[12 + k1]*dm4[120 + k2] + 2*dc[15 + k2]*dm4[123 + k1] + 2*dc[15 + k1]*dm4[123 + k2] + dc[24 + k2]*dm4[132 + k1] + dc[24 + k1]*dm4[132 + k2]);
      rhs(8) = -(a*ddm2[45 + 3*k1 + k2] + b[0]*ddm3[45 + 3*k1 + k2] + b[1]*ddm3[126 + 3*k1 + k2] + b[2]*ddm3[153 + 3*k1 + k2] + c[0]*ddm4[45 + 3*k1 + k2] + 2*c[1]*ddm4[126 + 3*k1 + k2] + 2*c[2]*ddm4[153 + 3*k1 + k2] + c[4]*ddm4[369 + 3*k1 + k2] + 2*c[5]*ddm4[396 + 3*k1 + k2] + c[8]*ddm4[477 + 3*k1 + k2] + da[k2]*dm2[15 + k1] + da[k1]*dm2[15 + k2] + db[k2]*dm3[15 + k1] + db[k1]*dm3[15 + k2] + db[3 + k2]*dm3[42 + k1] + db[3 + k1]*dm3[42 + k2] + db[6 + k2]*dm3[51 + k1] + db[6 + k1]*dm3[51 + k2] + dc[k2]*dm4[15 + k1] + dc[k1]*dm4[15 + k2] + 2*dc[3 + k2]*dm4[42 + k1] + 2*dc[3 + k1]*dm4[42 + k2] + 2*dc[6 + k2]*dm4[51 + k1] + 2*dc[6 + k1]*dm4[51 + k2] + dc[12 + k2]*dm4[123 + k1] + dc[12 + k1]*dm4[123 + k2] + 2*dc[15 + k2]*dm4[132 + k1] + 2*dc[15 + k1]*dm4[132 + k2] + dc[24 + k2]*dm4[159 + k1] + dc[24 + k1]*dm4[159 + k2]);
      rhs(9) = -(a*ddm2[72 + 3*k1 + k2] + b[0]*ddm3[72 + 3*k1 + k2] + b[1]*ddm3[153 + 3*k1 + k2] + b[2]*ddm3[234 + 3*k1 + k2] + c[0]*ddm4[72 + 3*k1 + k2] + 2*c[1]*ddm4[153 + 3*k1 + k2] + 2*c[2]*ddm4[234 + 3*k1 + k2] + c[4]*ddm4[396 + 3*k1 + k2] + 2*c[5]*ddm4[477 + 3*k1 + k2] + c[8]*ddm4[720 + 3*k1 + k2] + da[k2]*dm2[24 + k1] + da[k1]*dm2[24 + k2] + db[k2]*dm3[24 + k1] + db[k1]*dm3[24 + k2] + db[3 + k2]*dm3[51 + k1] + db[3 + k1]*dm3[51 + k2] + db[6 + k2]*dm3[78 + k1] + db[6 + k1]*dm3[78 + k2] + dc[k2]*dm4[24 + k1] + dc[k1]*dm4[24 + k2] + 2*dc[3 + k2]*dm4[51 + k1] + 2*dc[3 + k1]*dm4[51 + k2] + 2*dc[6 + k2]*dm4[78 + k1] + 2*dc[6 + k1]*dm4[78 + k2] + dc[12 + k2]*dm4[132 + k1] + dc[12 + k1]*dm4[132 + k2] + 2*dc[15 + k2]*dm4[159 + k1] + 2*dc[15 + k1]*dm4[159 + k2] + dc[24 + k2]*dm4[240 + k1] + dc[24 + k1]*dm4[240 + k2]);
 
      lhs = solver.solve(rhs);

      dda[3*k1 + k2] = lhs(0);
      ddb[3*k1 + k2] = lhs(1);
      ddb[9 + 3*k1 + k2] = lhs(2);
      ddb[18 + 3*k1 + k2] = lhs(3);
      ddc[3*k1 + k2] = lhs(4);
      ddc[9 + 3*k1 + k2] = lhs(5);
      ddc[18 + 3*k1 + k2] = lhs(6);
      ddc[36 + 3*k1 + k2] = lhs(7);
      ddc[45 + 3*k1 + k2] = lhs(8);
      ddc[72 + 3*k1 + k2] = lhs(9);
    }
  }
  
  // Fill symmetries
  fillSymmetries<Dim<3>>(ddb);
  fillSymmetries<Dim<3>>(c, dc, ddc);
}
// d=1, o=3
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<1>>& rkMoments,
                   Dim<1>::Scalar& a,
                   Dim<1>::Vector& b,
                   Dim<1>::Tensor& c,
                   Dim<1>::ThirdRankTensor& d,
                   Dim<1>::Vector& da,
                   Dim<1>::Tensor& db,
                   Dim<1>::ThirdRankTensor& dc,
                   Dim<1>::FourthRankTensor& dd,
                   Dim<1>::Tensor& dda,
                   Dim<1>::ThirdRankTensor& ddb,
                   Dim<1>::FourthRankTensor& ddc, 
                   Dim<1>::FifthRankTensor& ddd) {
  const auto dim = Dim<1>::nDim;
  const auto size = 4;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& m3 = rkMoments.m3;
  const auto& m4 = rkMoments.m4;
  const auto& m5 = rkMoments.m5;
  const auto& m6 = rkMoments.m6;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& dm3 = rkMoments.dm3;
  const auto& dm4 = rkMoments.dm4;
  const auto& dm5 = rkMoments.dm5;
  const auto& dm6 = rkMoments.dm6;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;
  const auto& ddm3 = rkMoments.ddm3;
  const auto& ddm4 = rkMoments.ddm4;
  const auto& ddm5 = rkMoments.ddm5;
  const auto& ddm6 = rkMoments.ddm6;
  const auto k1 = 0;
  const auto k2 = 0;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix << 
    m0, m1[0], m2[0], m3[0], 
    m1[0], m2[0], m3[0], m4[0], 
    m2[0], m3[0], m4[0], m5[0], 
    m3[0], m4[0], m5[0], m6[0];
  
  auto solver = matrix.colPivHouseholderQr();

  // Solve for values
  rhs(0) = 1;
  rhs(1) = 0;
  rhs(2) = 0;
  rhs(3) = 0;
  
  lhs = solver.solve(rhs);

  a = lhs(0);
  b[0] = lhs(1);
  c[0] = lhs(2);
  d[0] = lhs(3);
  
  // Solve for derivatives
  rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + c[0]*dm2[k1] + d[0]*dm3[k1]);
  rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + c[0]*dm3[k1] + d[0]*dm4[k1]);
  rhs(2) = -(a*dm2[k1] + b[0]*dm3[k1] + c[0]*dm4[k1] + d[0]*dm5[k1]);
  rhs(3) = -(a*dm3[k1] + b[0]*dm4[k1] + c[0]*dm5[k1] + d[0]*dm6[k1]);
  
  lhs = solver.solve(rhs);
  
  da[k1] = lhs(0);
  db[k1] = lhs(1);
  dc[k1] = lhs(2);
  dd[k1] = lhs(3);
 
  // Solve for second derivatives
  rhs(0) = -(a*ddm0[k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[k1 + k2] + c[0]*ddm2[k1 + k2] + d[0]*ddm3[k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + dc[k2]*dm2[k1] + dc[k1]*dm2[k2] + dd[k2]*dm3[k1] + dd[k1]*dm3[k2]);
  rhs(1) = -(a*ddm1[k1 + k2] + b[0]*ddm2[k1 + k2] + c[0]*ddm3[k1 + k2] + d[0]*ddm4[k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + dc[k2]*dm3[k1] + dc[k1]*dm3[k2] + dd[k2]*dm4[k1] + dd[k1]*dm4[k2]);
  rhs(2) = -(a*ddm2[k1 + k2] + b[0]*ddm3[k1 + k2] + c[0]*ddm4[k1 + k2] + d[0]*ddm5[k1 + k2] + da[k2]*dm2[k1] + da[k1]*dm2[k2] + db[k2]*dm3[k1] + db[k1]*dm3[k2] + dc[k2]*dm4[k1] + dc[k1]*dm4[k2] + dd[k2]*dm5[k1] + dd[k1]*dm5[k2]);
  rhs(3) = -(a*ddm3[k1 + k2] + b[0]*ddm4[k1 + k2] + c[0]*ddm5[k1 + k2] + d[0]*ddm6[k1 + k2] + da[k2]*dm3[k1] + da[k1]*dm3[k2] + db[k2]*dm4[k1] + db[k1]*dm4[k2] + dc[k2]*dm5[k1] + dc[k1]*dm5[k2] + dd[k2]*dm6[k1] + dd[k1]*dm6[k2]);

  lhs = solver.solve(rhs);

  dda[k1 + k2] = lhs(0);
  ddb[k1 + k2] = lhs(1);
  ddc[k1 + k2] = lhs(2);
  ddd[k1 + k2] = lhs(3);
}
// d=2, o=3
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<2>>& rkMoments,
                   Dim<2>::Scalar& a,
                   Dim<2>::Vector& b,
                   Dim<2>::Tensor& c,
                   Dim<2>::ThirdRankTensor& d,
                   Dim<2>::Vector& da,
                   Dim<2>::Tensor& db,
                   Dim<2>::ThirdRankTensor& dc,
                   Dim<2>::FourthRankTensor& dd,
                   Dim<2>::Tensor& dda,
                   Dim<2>::ThirdRankTensor& ddb,
                   Dim<2>::FourthRankTensor& ddc, 
                   Dim<2>::FifthRankTensor& ddd) {
  const auto dim = Dim<2>::nDim;
  const auto size = 10;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& m3 = rkMoments.m3;
  const auto& m4 = rkMoments.m4;
  const auto& m5 = rkMoments.m5;
  const auto& m6 = rkMoments.m6;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& dm3 = rkMoments.dm3;
  const auto& dm4 = rkMoments.dm4;
  const auto& dm5 = rkMoments.dm5;
  const auto& dm6 = rkMoments.dm6;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;
  const auto& ddm3 = rkMoments.ddm3;
  const auto& ddm4 = rkMoments.ddm4;
  const auto& ddm5 = rkMoments.ddm5;
  const auto& ddm6 = rkMoments.ddm6;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix <<
    m0, m1[0], m1[1], m2[0], 2*m2[1], m2[3], m3[0], 3*m3[1], 3*m3[3], m3[7], 
    m1[0], m2[0], m2[1], m3[0], 2*m3[1], m3[3], m4[0], 3*m4[1], 3*m4[3], m4[7], 
    m1[1], m2[1], m2[3], m3[1], 2*m3[3], m3[7], m4[1], 3*m4[3], 3*m4[7], m4[15], 
    m2[0], m3[0], m3[1], m4[0], 2*m4[1], m4[3], m5[0], 3*m5[1], 3*m5[3], m5[7], 
    m2[1], m3[1], m3[3], m4[1], 2*m4[3], m4[7], m5[1], 3*m5[3], 3*m5[7], m5[15], 
    m2[3], m3[3], m3[7], m4[3], 2*m4[7], m4[15], m5[3], 3*m5[7], 3*m5[15], m5[31], 
    m3[0], m4[0], m4[1], m5[0], 2*m5[1], m5[3], m6[0], 3*m6[1], 3*m6[3], m6[7], 
    m3[1], m4[1], m4[3], m5[1], 2*m5[3], m5[7], m6[1], 3*m6[3], 3*m6[7], m6[15], 
    m3[3], m4[3], m4[7], m5[3], 2*m5[7], m5[15], m6[3], 3*m6[7], 3*m6[15], m6[31], 
    m3[7], m4[7], m4[15], m5[7], 2*m5[15], m5[31], m6[7], 3*m6[15], 3*m6[31], m6[63];
  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = 1;
  rhs(1) = 0;
  rhs(2) = 0;
  rhs(3) = 0;
  rhs(4) = 0;
  rhs(5) = 0;
  rhs(6) = 0;
  rhs(7) = 0;
  rhs(8) = 0;
  rhs(9) = 0;
 
  lhs = solver.solve(rhs);

  a = lhs(0);
  b[0] = lhs(1);
  b[1] = lhs(2);
  c[0] = lhs(3);
  c[1] = lhs(4);
  c[3] = lhs(5);
  d[0] = lhs(6);
  d[1] = lhs(7);
  d[3] = lhs(8);
  d[7] = lhs(9);
 
  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[2 + k1] + c[0]*dm2[k1] + 2*c[1]*dm2[2 + k1] + c[3]*dm2[6 + k1] + d[0]*dm3[k1] + 3*d[1]*dm3[2 + k1] + 3*d[3]*dm3[6 + k1] + d[7]*dm3[14 + k1]);
    rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + b[1]*dm2[2 + k1] + c[0]*dm3[k1] + 2*c[1]*dm3[2 + k1] + c[3]*dm3[6 + k1] + d[0]*dm4[k1] + 3*d[1]*dm4[2 + k1] + 3*d[3]*dm4[6 + k1] + d[7]*dm4[14 + k1]);
    rhs(2) = -(a*dm1[2 + k1] + b[0]*dm2[2 + k1] + b[1]*dm2[6 + k1] + c[0]*dm3[2 + k1] + 2*c[1]*dm3[6 + k1] + c[3]*dm3[14 + k1] + d[0]*dm4[2 + k1] + 3*d[1]*dm4[6 + k1] + 3*d[3]*dm4[14 + k1] + d[7]*dm4[30 + k1]);
    rhs(3) = -(a*dm2[k1] + b[0]*dm3[k1] + b[1]*dm3[2 + k1] + c[0]*dm4[k1] + 2*c[1]*dm4[2 + k1] + c[3]*dm4[6 + k1] + d[0]*dm5[k1] + 3*d[1]*dm5[2 + k1] + 3*d[3]*dm5[6 + k1] + d[7]*dm5[14 + k1]);
    rhs(4) = -(a*dm2[2 + k1] + b[0]*dm3[2 + k1] + b[1]*dm3[6 + k1] + c[0]*dm4[2 + k1] + 2*c[1]*dm4[6 + k1] + c[3]*dm4[14 + k1] + d[0]*dm5[2 + k1] + 3*d[1]*dm5[6 + k1] + 3*d[3]*dm5[14 + k1] + d[7]*dm5[30 + k1]);
    rhs(5) = -(a*dm2[6 + k1] + b[0]*dm3[6 + k1] + b[1]*dm3[14 + k1] + c[0]*dm4[6 + k1] + 2*c[1]*dm4[14 + k1] + c[3]*dm4[30 + k1] + d[0]*dm5[6 + k1] + 3*d[1]*dm5[14 + k1] + 3*d[3]*dm5[30 + k1] + d[7]*dm5[62 + k1]);
    rhs(6) = -(a*dm3[k1] + b[0]*dm4[k1] + b[1]*dm4[2 + k1] + c[0]*dm5[k1] + 2*c[1]*dm5[2 + k1] + c[3]*dm5[6 + k1] + d[0]*dm6[k1] + 3*d[1]*dm6[2 + k1] + 3*d[3]*dm6[6 + k1] + d[7]*dm6[14 + k1]);
    rhs(7) = -(a*dm3[2 + k1] + b[0]*dm4[2 + k1] + b[1]*dm4[6 + k1] + c[0]*dm5[2 + k1] + 2*c[1]*dm5[6 + k1] + c[3]*dm5[14 + k1] + d[0]*dm6[2 + k1] + 3*d[1]*dm6[6 + k1] + 3*d[3]*dm6[14 + k1] + d[7]*dm6[30 + k1]);
    rhs(8) = -(a*dm3[6 + k1] + b[0]*dm4[6 + k1] + b[1]*dm4[14 + k1] + c[0]*dm5[6 + k1] + 2*c[1]*dm5[14 + k1] + c[3]*dm5[30 + k1] + d[0]*dm6[6 + k1] + 3*d[1]*dm6[14 + k1] + 3*d[3]*dm6[30 + k1] + d[7]*dm6[62 + k1]);
    rhs(9) = -(a*dm3[14 + k1] + b[0]*dm4[14 + k1] + b[1]*dm4[30 + k1] + c[0]*dm5[14 + k1] + 2*c[1]*dm5[30 + k1] + c[3]*dm5[62 + k1] + d[0]*dm6[14 + k1] + 3*d[1]*dm6[30 + k1] + 3*d[3]*dm6[62 + k1] + d[7]*dm6[126 + k1]);
    
    lhs = solver.solve(rhs);

    da[k1] = lhs(0);
    db[k1] = lhs(1);
    db[2 + k1] = lhs(2);
    dc[k1] = lhs(3);
    dc[2 + k1] = lhs(4);
    dc[6 + k1] = lhs(5);
    dd[k1] = lhs(6);
    dd[2 + k1] = lhs(7);
    dd[6 + k1] = lhs(8);
    dd[14 + k1] = lhs(9);    
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[2*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[2*k1 + k2] + b[1]*ddm1[4 + 2*k1 + k2] + c[0]*ddm2[2*k1 + k2] + 2*c[1]*ddm2[4 + 2*k1 + k2] + c[3]*ddm2[12 + 2*k1 + k2] + d[0]*ddm3[2*k1 + k2] + 3*d[1]*ddm3[4 + 2*k1 + k2] + 3*d[3]*ddm3[12 + 2*k1 + k2] + d[7]*ddm3[28 + 2*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[2 + k2]*dm1[2 + k1] + db[2 + k1]*dm1[2 + k2] + dc[k2]*dm2[k1] + dc[k1]*dm2[k2] + 2*dc[2 + k2]*dm2[2 + k1] + 2*dc[2 + k1]*dm2[2 + k2] + dc[6 + k2]*dm2[6 + k1] + dc[6 + k1]*dm2[6 + k2] + dd[k2]*dm3[k1] + dd[k1]*dm3[k2] + 3*dd[2 + k2]*dm3[2 + k1] + 3*dd[2 + k1]*dm3[2 + k2] + 3*dd[6 + k2]*dm3[6 + k1] + 3*dd[6 + k1]*dm3[6 + k2] + dd[14 + k2]*dm3[14 + k1] + dd[14 + k1]*dm3[14 + k2]);
      rhs(1) = -(a*ddm1[2*k1 + k2] + b[0]*ddm2[2*k1 + k2] + b[1]*ddm2[4 + 2*k1 + k2] + c[0]*ddm3[2*k1 + k2] + 2*c[1]*ddm3[4 + 2*k1 + k2] + c[3]*ddm3[12 + 2*k1 + k2] + d[0]*ddm4[2*k1 + k2] + 3*d[1]*ddm4[4 + 2*k1 + k2] + 3*d[3]*ddm4[12 + 2*k1 + k2] + d[7]*ddm4[28 + 2*k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + db[2 + k2]*dm2[2 + k1] + db[2 + k1]*dm2[2 + k2] + dc[k2]*dm3[k1] + dc[k1]*dm3[k2] + 2*dc[2 + k2]*dm3[2 + k1] + 2*dc[2 + k1]*dm3[2 + k2] + dc[6 + k2]*dm3[6 + k1] + dc[6 + k1]*dm3[6 + k2] + dd[k2]*dm4[k1] + dd[k1]*dm4[k2] + 3*dd[2 + k2]*dm4[2 + k1] + 3*dd[2 + k1]*dm4[2 + k2] + 3*dd[6 + k2]*dm4[6 + k1] + 3*dd[6 + k1]*dm4[6 + k2] + dd[14 + k2]*dm4[14 + k1] + dd[14 + k1]*dm4[14 + k2]);
      rhs(2) = -(a*ddm1[4 + 2*k1 + k2] + b[0]*ddm2[4 + 2*k1 + k2] + b[1]*ddm2[12 + 2*k1 + k2] + c[0]*ddm3[4 + 2*k1 + k2] + 2*c[1]*ddm3[12 + 2*k1 + k2] + c[3]*ddm3[28 + 2*k1 + k2] + d[0]*ddm4[4 + 2*k1 + k2] + 3*d[1]*ddm4[12 + 2*k1 + k2] + 3*d[3]*ddm4[28 + 2*k1 + k2] + d[7]*ddm4[60 + 2*k1 + k2] + da[k2]*dm1[2 + k1] + da[k1]*dm1[2 + k2] + db[k2]*dm2[2 + k1] + db[k1]*dm2[2 + k2] + db[2 + k2]*dm2[6 + k1] + db[2 + k1]*dm2[6 + k2] + dc[k2]*dm3[2 + k1] + dc[k1]*dm3[2 + k2] + 2*dc[2 + k2]*dm3[6 + k1] + 2*dc[2 + k1]*dm3[6 + k2] + dc[6 + k2]*dm3[14 + k1] + dc[6 + k1]*dm3[14 + k2] + dd[k2]*dm4[2 + k1] + dd[k1]*dm4[2 + k2] + 3*dd[2 + k2]*dm4[6 + k1] + 3*dd[2 + k1]*dm4[6 + k2] + 3*dd[6 + k2]*dm4[14 + k1] + 3*dd[6 + k1]*dm4[14 + k2] + dd[14 + k2]*dm4[30 + k1] + dd[14 + k1]*dm4[30 + k2]);
      rhs(3) = -(a*ddm2[2*k1 + k2] + b[0]*ddm3[2*k1 + k2] + b[1]*ddm3[4 + 2*k1 + k2] + c[0]*ddm4[2*k1 + k2] + 2*c[1]*ddm4[4 + 2*k1 + k2] + c[3]*ddm4[12 + 2*k1 + k2] + d[0]*ddm5[2*k1 + k2] + 3*d[1]*ddm5[4 + 2*k1 + k2] + 3*d[3]*ddm5[12 + 2*k1 + k2] + d[7]*ddm5[28 + 2*k1 + k2] + da[k2]*dm2[k1] + da[k1]*dm2[k2] + db[k2]*dm3[k1] + db[k1]*dm3[k2] + db[2 + k2]*dm3[2 + k1] + db[2 + k1]*dm3[2 + k2] + dc[k2]*dm4[k1] + dc[k1]*dm4[k2] + 2*dc[2 + k2]*dm4[2 + k1] + 2*dc[2 + k1]*dm4[2 + k2] + dc[6 + k2]*dm4[6 + k1] + dc[6 + k1]*dm4[6 + k2] + dd[k2]*dm5[k1] + dd[k1]*dm5[k2] + 3*dd[2 + k2]*dm5[2 + k1] + 3*dd[2 + k1]*dm5[2 + k2] + 3*dd[6 + k2]*dm5[6 + k1] + 3*dd[6 + k1]*dm5[6 + k2] + dd[14 + k2]*dm5[14 + k1] + dd[14 + k1]*dm5[14 + k2]);
      rhs(4) = -(a*ddm2[4 + 2*k1 + k2] + b[0]*ddm3[4 + 2*k1 + k2] + b[1]*ddm3[12 + 2*k1 + k2] + c[0]*ddm4[4 + 2*k1 + k2] + 2*c[1]*ddm4[12 + 2*k1 + k2] + c[3]*ddm4[28 + 2*k1 + k2] + d[0]*ddm5[4 + 2*k1 + k2] + 3*d[1]*ddm5[12 + 2*k1 + k2] + 3*d[3]*ddm5[28 + 2*k1 + k2] + d[7]*ddm5[60 + 2*k1 + k2] + da[k2]*dm2[2 + k1] + da[k1]*dm2[2 + k2] + db[k2]*dm3[2 + k1] + db[k1]*dm3[2 + k2] + db[2 + k2]*dm3[6 + k1] + db[2 + k1]*dm3[6 + k2] + dc[k2]*dm4[2 + k1] + dc[k1]*dm4[2 + k2] + 2*dc[2 + k2]*dm4[6 + k1] + 2*dc[2 + k1]*dm4[6 + k2] + dc[6 + k2]*dm4[14 + k1] + dc[6 + k1]*dm4[14 + k2] + dd[k2]*dm5[2 + k1] + dd[k1]*dm5[2 + k2] + 3*dd[2 + k2]*dm5[6 + k1] + 3*dd[2 + k1]*dm5[6 + k2] + 3*dd[6 + k2]*dm5[14 + k1] + 3*dd[6 + k1]*dm5[14 + k2] + dd[14 + k2]*dm5[30 + k1] + dd[14 + k1]*dm5[30 + k2]);
      rhs(5) = -(a*ddm2[12 + 2*k1 + k2] + b[0]*ddm3[12 + 2*k1 + k2] + b[1]*ddm3[28 + 2*k1 + k2] + c[0]*ddm4[12 + 2*k1 + k2] + 2*c[1]*ddm4[28 + 2*k1 + k2] + c[3]*ddm4[60 + 2*k1 + k2] + d[0]*ddm5[12 + 2*k1 + k2] + 3*d[1]*ddm5[28 + 2*k1 + k2] + 3*d[3]*ddm5[60 + 2*k1 + k2] + d[7]*ddm5[124 + 2*k1 + k2] + da[k2]*dm2[6 + k1] + da[k1]*dm2[6 + k2] + db[k2]*dm3[6 + k1] + db[k1]*dm3[6 + k2] + db[2 + k2]*dm3[14 + k1] + db[2 + k1]*dm3[14 + k2] + dc[k2]*dm4[6 + k1] + dc[k1]*dm4[6 + k2] + 2*dc[2 + k2]*dm4[14 + k1] + 2*dc[2 + k1]*dm4[14 + k2] + dc[6 + k2]*dm4[30 + k1] + dc[6 + k1]*dm4[30 + k2] + dd[k2]*dm5[6 + k1] + dd[k1]*dm5[6 + k2] + 3*dd[2 + k2]*dm5[14 + k1] + 3*dd[2 + k1]*dm5[14 + k2] + 3*dd[6 + k2]*dm5[30 + k1] + 3*dd[6 + k1]*dm5[30 + k2] + dd[14 + k2]*dm5[62 + k1] + dd[14 + k1]*dm5[62 + k2]);
      rhs(6) = -(a*ddm3[2*k1 + k2] + b[0]*ddm4[2*k1 + k2] + b[1]*ddm4[4 + 2*k1 + k2] + c[0]*ddm5[2*k1 + k2] + 2*c[1]*ddm5[4 + 2*k1 + k2] + c[3]*ddm5[12 + 2*k1 + k2] + d[0]*ddm6[2*k1 + k2] + 3*d[1]*ddm6[4 + 2*k1 + k2] + 3*d[3]*ddm6[12 + 2*k1 + k2] + d[7]*ddm6[28 + 2*k1 + k2] + da[k2]*dm3[k1] + da[k1]*dm3[k2] + db[k2]*dm4[k1] + db[k1]*dm4[k2] + db[2 + k2]*dm4[2 + k1] + db[2 + k1]*dm4[2 + k2] + dc[k2]*dm5[k1] + dc[k1]*dm5[k2] + 2*dc[2 + k2]*dm5[2 + k1] + 2*dc[2 + k1]*dm5[2 + k2] + dc[6 + k2]*dm5[6 + k1] + dc[6 + k1]*dm5[6 + k2] + dd[k2]*dm6[k1] + dd[k1]*dm6[k2] + 3*dd[2 + k2]*dm6[2 + k1] + 3*dd[2 + k1]*dm6[2 + k2] + 3*dd[6 + k2]*dm6[6 + k1] + 3*dd[6 + k1]*dm6[6 + k2] + dd[14 + k2]*dm6[14 + k1] + dd[14 + k1]*dm6[14 + k2]);
      rhs(7) = -(a*ddm3[4 + 2*k1 + k2] + b[0]*ddm4[4 + 2*k1 + k2] + b[1]*ddm4[12 + 2*k1 + k2] + c[0]*ddm5[4 + 2*k1 + k2] + 2*c[1]*ddm5[12 + 2*k1 + k2] + c[3]*ddm5[28 + 2*k1 + k2] + d[0]*ddm6[4 + 2*k1 + k2] + 3*d[1]*ddm6[12 + 2*k1 + k2] + 3*d[3]*ddm6[28 + 2*k1 + k2] + d[7]*ddm6[60 + 2*k1 + k2] + da[k2]*dm3[2 + k1] + da[k1]*dm3[2 + k2] + db[k2]*dm4[2 + k1] + db[k1]*dm4[2 + k2] + db[2 + k2]*dm4[6 + k1] + db[2 + k1]*dm4[6 + k2] + dc[k2]*dm5[2 + k1] + dc[k1]*dm5[2 + k2] + 2*dc[2 + k2]*dm5[6 + k1] + 2*dc[2 + k1]*dm5[6 + k2] + dc[6 + k2]*dm5[14 + k1] + dc[6 + k1]*dm5[14 + k2] + dd[k2]*dm6[2 + k1] + dd[k1]*dm6[2 + k2] + 3*dd[2 + k2]*dm6[6 + k1] + 3*dd[2 + k1]*dm6[6 + k2] + 3*dd[6 + k2]*dm6[14 + k1] + 3*dd[6 + k1]*dm6[14 + k2] + dd[14 + k2]*dm6[30 + k1] + dd[14 + k1]*dm6[30 + k2]);
      rhs(8) = -(a*ddm3[12 + 2*k1 + k2] + b[0]*ddm4[12 + 2*k1 + k2] + b[1]*ddm4[28 + 2*k1 + k2] + c[0]*ddm5[12 + 2*k1 + k2] + 2*c[1]*ddm5[28 + 2*k1 + k2] + c[3]*ddm5[60 + 2*k1 + k2] + d[0]*ddm6[12 + 2*k1 + k2] + 3*d[1]*ddm6[28 + 2*k1 + k2] + 3*d[3]*ddm6[60 + 2*k1 + k2] + d[7]*ddm6[124 + 2*k1 + k2] + da[k2]*dm3[6 + k1] + da[k1]*dm3[6 + k2] + db[k2]*dm4[6 + k1] + db[k1]*dm4[6 + k2] + db[2 + k2]*dm4[14 + k1] + db[2 + k1]*dm4[14 + k2] + dc[k2]*dm5[6 + k1] + dc[k1]*dm5[6 + k2] + 2*dc[2 + k2]*dm5[14 + k1] + 2*dc[2 + k1]*dm5[14 + k2] + dc[6 + k2]*dm5[30 + k1] + dc[6 + k1]*dm5[30 + k2] + dd[k2]*dm6[6 + k1] + dd[k1]*dm6[6 + k2] + 3*dd[2 + k2]*dm6[14 + k1] + 3*dd[2 + k1]*dm6[14 + k2] + 3*dd[6 + k2]*dm6[30 + k1] + 3*dd[6 + k1]*dm6[30 + k2] + dd[14 + k2]*dm6[62 + k1] + dd[14 + k1]*dm6[62 + k2]);
      rhs(9) = -(a*ddm3[28 + 2*k1 + k2] + b[0]*ddm4[28 + 2*k1 + k2] + b[1]*ddm4[60 + 2*k1 + k2] + c[0]*ddm5[28 + 2*k1 + k2] + 2*c[1]*ddm5[60 + 2*k1 + k2] + c[3]*ddm5[124 + 2*k1 + k2] + d[0]*ddm6[28 + 2*k1 + k2] + 3*d[1]*ddm6[60 + 2*k1 + k2] + 3*d[3]*ddm6[124 + 2*k1 + k2] + d[7]*ddm6[252 + 2*k1 + k2] + da[k2]*dm3[14 + k1] + da[k1]*dm3[14 + k2] + db[k2]*dm4[14 + k1] + db[k1]*dm4[14 + k2] + db[2 + k2]*dm4[30 + k1] + db[2 + k1]*dm4[30 + k2] + dc[k2]*dm5[14 + k1] + dc[k1]*dm5[14 + k2] + 2*dc[2 + k2]*dm5[30 + k1] + 2*dc[2 + k1]*dm5[30 + k2] + dc[6 + k2]*dm5[62 + k1] + dc[6 + k1]*dm5[62 + k2] + dd[k2]*dm6[14 + k1] + dd[k1]*dm6[14 + k2] + 3*dd[2 + k2]*dm6[30 + k1] + 3*dd[2 + k1]*dm6[30 + k2] + 3*dd[6 + k2]*dm6[62 + k1] + 3*dd[6 + k1]*dm6[62 + k2] + dd[14 + k2]*dm6[126 + k1] + dd[14 + k1]*dm6[126 + k2]);
 
      lhs = solver.solve(rhs);

      dda[2*k1 + k2] = lhs(0);
      ddb[2*k1 + k2] = lhs(1);
      ddb[4 + 2*k1 + k2] = lhs(2);
      ddc[2*k1 + k2] = lhs(3);
      ddc[4 + 2*k1 + k2] = lhs(4);
      ddc[12 + 2*k1 + k2] = lhs(5);
      ddd[2*k1 + k2] = lhs(6);
      ddd[4 + 2*k1 + k2] = lhs(7);
      ddd[12 + 2*k1 + k2] = lhs(8);
      ddd[28 + 2*k1 + k2] = lhs(9);      
    }
  }
  
  // Fill symmetries
  fillSymmetries<Dim<2>>(ddb);
  fillSymmetries<Dim<2>>(c, dc, ddc);
  fillSymmetries<Dim<2>>(d, dd, ddd);
}
// d=3, o=3
template<>
inline
void
computeCorrections(const RKMomentValues<Dim<3>>& rkMoments,
                   Dim<3>::Scalar& a,
                   Dim<3>::Vector& b,
                   Dim<3>::Tensor& c,
                   Dim<3>::ThirdRankTensor& d,
                   Dim<3>::Vector& da,
                   Dim<3>::Tensor& db,
                   Dim<3>::ThirdRankTensor& dc,
                   Dim<3>::FourthRankTensor& dd,
                   Dim<3>::Tensor& dda,
                   Dim<3>::ThirdRankTensor& ddb,
                   Dim<3>::FourthRankTensor& ddc, 
                   Dim<3>::FifthRankTensor& ddd) {
  const auto dim = Dim<3>::nDim;
  const auto size = 20;
  typedef Eigen::Matrix<double, size, size> MatrixType;
  typedef Eigen::Matrix<double, size, 1> VectorType;
  const auto& m0 = rkMoments.m0;
  const auto& m1 = rkMoments.m1;
  const auto& m2 = rkMoments.m2;
  const auto& m3 = rkMoments.m3;
  const auto& m4 = rkMoments.m4;
  const auto& m5 = rkMoments.m5;
  const auto& m6 = rkMoments.m6;
  const auto& dm0 = rkMoments.dm0;
  const auto& dm1 = rkMoments.dm1;
  const auto& dm2 = rkMoments.dm2;
  const auto& dm3 = rkMoments.dm3;
  const auto& dm4 = rkMoments.dm4;
  const auto& dm5 = rkMoments.dm5;
  const auto& dm6 = rkMoments.dm6;
  const auto& ddm0 = rkMoments.ddm0;
  const auto& ddm1 = rkMoments.ddm1;
  const auto& ddm2 = rkMoments.ddm2;
  const auto& ddm3 = rkMoments.ddm3;
  const auto& ddm4 = rkMoments.ddm4;
  const auto& ddm5 = rkMoments.ddm5;
  const auto& ddm6 = rkMoments.ddm6;

  // Get matrix to invert
  VectorType rhs;
  VectorType lhs;
  MatrixType matrix;
  matrix <<
    m0, m1[0], m1[1], m1[2], m2[0], 2*m2[1], 2*m2[2], m2[4], 2*m2[5], m2[8], m3[0], 3*m3[1], 3*m3[2], 3*m3[4], 6*m3[5], 3*m3[8], m3[13], 3*m3[14], 3*m3[17], m3[26], 
    m1[0], m2[0], m2[1], m2[2], m3[0], 2*m3[1], 2*m3[2], m3[4], 2*m3[5], m3[8], m4[0], 3*m4[1], 3*m4[2], 3*m4[4], 6*m4[5], 3*m4[8], m4[13], 3*m4[14], 3*m4[17], m4[26], 
    m1[1], m2[1], m2[4], m2[5], m3[1], 2*m3[4], 2*m3[5], m3[13], 2*m3[14], m3[17], m4[1], 3*m4[4], 3*m4[5], 3*m4[13], 6*m4[14], 3*m4[17], m4[40], 3*m4[41], 3*m4[44], m4[53], 
    m1[2], m2[2], m2[5], m2[8], m3[2], 2*m3[5], 2*m3[8], m3[14], 2*m3[17], m3[26], m4[2], 3*m4[5], 3*m4[8], 3*m4[14], 6*m4[17], 3*m4[26], m4[41], 3*m4[44], 3*m4[53], m4[80], 
    m2[0], m3[0], m3[1], m3[2], m4[0], 2*m4[1], 2*m4[2], m4[4], 2*m4[5], m4[8], m5[0], 3*m5[1], 3*m5[2], 3*m5[4], 6*m5[5], 3*m5[8], m5[13], 3*m5[14], 3*m5[17], m5[26], 
    m2[1], m3[1], m3[4], m3[5], m4[1], 2*m4[4], 2*m4[5], m4[13], 2*m4[14], m4[17], m5[1], 3*m5[4], 3*m5[5], 3*m5[13], 6*m5[14], 3*m5[17], m5[40], 3*m5[41], 3*m5[44], m5[53], 
    m2[2], m3[2], m3[5], m3[8], m4[2], 2*m4[5], 2*m4[8], m4[14], 2*m4[17], m4[26], m5[2], 3*m5[5], 3*m5[8], 3*m5[14], 6*m5[17], 3*m5[26], m5[41], 3*m5[44], 3*m5[53], m5[80], 
    m2[4], m3[4], m3[13], m3[14], m4[4], 2*m4[13], 2*m4[14], m4[40], 2*m4[41], m4[44], m5[4], 3*m5[13], 3*m5[14], 3*m5[40], 6*m5[41], 3*m5[44], m5[121], 3*m5[122], 3*m5[125], m5[134], 
    m2[5], m3[5], m3[14], m3[17], m4[5], 2*m4[14], 2*m4[17], m4[41], 2*m4[44], m4[53], m5[5], 3*m5[14], 3*m5[17], 3*m5[41], 6*m5[44], 3*m5[53], m5[122], 3*m5[125], 3*m5[134], m5[161], 
    m2[8], m3[8], m3[17], m3[26], m4[8], 2*m4[17], 2*m4[26], m4[44], 2*m4[53], m4[80], m5[8], 3*m5[17], 3*m5[26], 3*m5[44], 6*m5[53], 3*m5[80], m5[125], 3*m5[134], 3*m5[161], m5[242], 
    m3[0], m4[0], m4[1], m4[2], m5[0], 2*m5[1], 2*m5[2], m5[4], 2*m5[5], m5[8], m6[0], 3*m6[1], 3*m6[2], 3*m6[4], 6*m6[5], 3*m6[8], m6[13], 3*m6[14], 3*m6[17], m6[26], 
    m3[1], m4[1], m4[4], m4[5], m5[1], 2*m5[4], 2*m5[5], m5[13], 2*m5[14], m5[17], m6[1], 3*m6[4], 3*m6[5], 3*m6[13], 6*m6[14], 3*m6[17], m6[40], 3*m6[41], 3*m6[44], m6[53], 
    m3[2], m4[2], m4[5], m4[8], m5[2], 2*m5[5], 2*m5[8], m5[14], 2*m5[17], m5[26], m6[2], 3*m6[5], 3*m6[8], 3*m6[14], 6*m6[17], 3*m6[26], m6[41], 3*m6[44], 3*m6[53], m6[80], 
    m3[4], m4[4], m4[13], m4[14], m5[4], 2*m5[13], 2*m5[14], m5[40], 2*m5[41], m5[44], m6[4], 3*m6[13], 3*m6[14], 3*m6[40], 6*m6[41], 3*m6[44], m6[121], 3*m6[122], 3*m6[125], m6[134], 
    m3[5], m4[5], m4[14], m4[17], m5[5], 2*m5[14], 2*m5[17], m5[41], 2*m5[44], m5[53], m6[5], 3*m6[14], 3*m6[17], 3*m6[41], 6*m6[44], 3*m6[53], m6[122], 3*m6[125], 3*m6[134], m6[161], 
    m3[8], m4[8], m4[17], m4[26], m5[8], 2*m5[17], 2*m5[26], m5[44], 2*m5[53], m5[80], m6[8], 3*m6[17], 3*m6[26], 3*m6[44], 6*m6[53], 3*m6[80], m6[125], 3*m6[134], 3*m6[161], m6[242], 
    m3[13], m4[13], m4[40], m4[41], m5[13], 2*m5[40], 2*m5[41], m5[121], 2*m5[122], m5[125], m6[13], 3*m6[40], 3*m6[41], 3*m6[121], 6*m6[122], 3*m6[125], m6[364], 3*m6[365], 3*m6[368], m6[377], 
    m3[14], m4[14], m4[41], m4[44], m5[14], 2*m5[41], 2*m5[44], m5[122], 2*m5[125], m5[134], m6[14], 3*m6[41], 3*m6[44], 3*m6[122], 6*m6[125], 3*m6[134], m6[365], 3*m6[368], 3*m6[377], m6[404], 
    m3[17], m4[17], m4[44], m4[53], m5[17], 2*m5[44], 2*m5[53], m5[125], 2*m5[134], m5[161], m6[17], 3*m6[44], 3*m6[53], 3*m6[125], 6*m6[134], 3*m6[161], m6[368], 3*m6[377], 3*m6[404], m6[485], 
    m3[26], m4[26], m4[53], m4[80], m5[26], 2*m5[53], 2*m5[80], m5[134], 2*m5[161], m5[242], m6[26], 3*m6[53], 3*m6[80], 3*m6[134], 6*m6[161], 3*m6[242], m6[377], 3*m6[404], 3*m6[485], m6[728];    
  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = 1;
  rhs(1) = 0;
  rhs(2) = 0;
  rhs(3) = 0;
  rhs(4) = 0;
  rhs(5) = 0;
  rhs(6) = 0;
  rhs(7) = 0;
  rhs(8) = 0;
  rhs(9) = 0;
  rhs(10) = 0;
  rhs(11) = 0;
  rhs(12) = 0;
  rhs(13) = 0;
  rhs(14) = 0;
  rhs(15) = 0;
  rhs(16) = 0;
  rhs(17) = 0;
  rhs(18) = 0;
  rhs(19) = 0;
  
  lhs = solver.solve(rhs);

  a = lhs(0);
  b[0] = lhs(1);
  b[1] = lhs(2);
  b[2] = lhs(3);
  c[0] = lhs(4);
  c[1] = lhs(5);
  c[2] = lhs(6);
  c[4] = lhs(7);
  c[5] = lhs(8);
  c[8] = lhs(9);
  d[0] = lhs(10);
  d[1] = lhs(11);
  d[2] = lhs(12);
  d[4] = lhs(13);
  d[5] = lhs(14);
  d[8] = lhs(15);
  d[13] = lhs(16);
  d[14] = lhs(17);
  d[17] = lhs(18);
  d[26] = lhs(19);

  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[3 + k1] + b[2]*dm1[6 + k1] + c[0]*dm2[k1] + 2*c[1]*dm2[3 + k1] + 2*c[2]*dm2[6 + k1] + c[4]*dm2[12 + k1] + 2*c[5]*dm2[15 + k1] + c[8]*dm2[24 + k1] + d[0]*dm3[k1] + 3*d[1]*dm3[3 + k1] + 3*d[2]*dm3[6 + k1] + 3*d[4]*dm3[12 + k1] + 6*d[5]*dm3[15 + k1] + 3*d[8]*dm3[24 + k1] + d[13]*dm3[39 + k1] + 3*d[14]*dm3[42 + k1] + 3*d[17]*dm3[51 + k1] + d[26]*dm3[78 + k1]);
    rhs(1) = -(a*dm1[k1] + b[0]*dm2[k1] + b[1]*dm2[3 + k1] + b[2]*dm2[6 + k1] + c[0]*dm3[k1] + 2*c[1]*dm3[3 + k1] + 2*c[2]*dm3[6 + k1] + c[4]*dm3[12 + k1] + 2*c[5]*dm3[15 + k1] + c[8]*dm3[24 + k1] + d[0]*dm4[k1] + 3*d[1]*dm4[3 + k1] + 3*d[2]*dm4[6 + k1] + 3*d[4]*dm4[12 + k1] + 6*d[5]*dm4[15 + k1] + 3*d[8]*dm4[24 + k1] + d[13]*dm4[39 + k1] + 3*d[14]*dm4[42 + k1] + 3*d[17]*dm4[51 + k1] + d[26]*dm4[78 + k1]);
    rhs(2) = -(a*dm1[3 + k1] + b[0]*dm2[3 + k1] + b[1]*dm2[12 + k1] + b[2]*dm2[15 + k1] + c[0]*dm3[3 + k1] + 2*c[1]*dm3[12 + k1] + 2*c[2]*dm3[15 + k1] + c[4]*dm3[39 + k1] + 2*c[5]*dm3[42 + k1] + c[8]*dm3[51 + k1] + d[0]*dm4[3 + k1] + 3*d[1]*dm4[12 + k1] + 3*d[2]*dm4[15 + k1] + 3*d[4]*dm4[39 + k1] + 6*d[5]*dm4[42 + k1] + 3*d[8]*dm4[51 + k1] + d[13]*dm4[120 + k1] + 3*d[14]*dm4[123 + k1] + 3*d[17]*dm4[132 + k1] + d[26]*dm4[159 + k1]);
    rhs(3) = -(a*dm1[6 + k1] + b[0]*dm2[6 + k1] + b[1]*dm2[15 + k1] + b[2]*dm2[24 + k1] + c[0]*dm3[6 + k1] + 2*c[1]*dm3[15 + k1] + 2*c[2]*dm3[24 + k1] + c[4]*dm3[42 + k1] + 2*c[5]*dm3[51 + k1] + c[8]*dm3[78 + k1] + d[0]*dm4[6 + k1] + 3*d[1]*dm4[15 + k1] + 3*d[2]*dm4[24 + k1] + 3*d[4]*dm4[42 + k1] + 6*d[5]*dm4[51 + k1] + 3*d[8]*dm4[78 + k1] + d[13]*dm4[123 + k1] + 3*d[14]*dm4[132 + k1] + 3*d[17]*dm4[159 + k1] + d[26]*dm4[240 + k1]);
    rhs(4) = -(a*dm2[k1] + b[0]*dm3[k1] + b[1]*dm3[3 + k1] + b[2]*dm3[6 + k1] + c[0]*dm4[k1] + 2*c[1]*dm4[3 + k1] + 2*c[2]*dm4[6 + k1] + c[4]*dm4[12 + k1] + 2*c[5]*dm4[15 + k1] + c[8]*dm4[24 + k1] + d[0]*dm5[k1] + 3*d[1]*dm5[3 + k1] + 3*d[2]*dm5[6 + k1] + 3*d[4]*dm5[12 + k1] + 6*d[5]*dm5[15 + k1] + 3*d[8]*dm5[24 + k1] + d[13]*dm5[39 + k1] + 3*d[14]*dm5[42 + k1] + 3*d[17]*dm5[51 + k1] + d[26]*dm5[78 + k1]);
    rhs(5) = -(a*dm2[3 + k1] + b[0]*dm3[3 + k1] + b[1]*dm3[12 + k1] + b[2]*dm3[15 + k1] + c[0]*dm4[3 + k1] + 2*c[1]*dm4[12 + k1] + 2*c[2]*dm4[15 + k1] + c[4]*dm4[39 + k1] + 2*c[5]*dm4[42 + k1] + c[8]*dm4[51 + k1] + d[0]*dm5[3 + k1] + 3*d[1]*dm5[12 + k1] + 3*d[2]*dm5[15 + k1] + 3*d[4]*dm5[39 + k1] + 6*d[5]*dm5[42 + k1] + 3*d[8]*dm5[51 + k1] + d[13]*dm5[120 + k1] + 3*d[14]*dm5[123 + k1] + 3*d[17]*dm5[132 + k1] + d[26]*dm5[159 + k1]);
    rhs(6) = -(a*dm2[6 + k1] + b[0]*dm3[6 + k1] + b[1]*dm3[15 + k1] + b[2]*dm3[24 + k1] + c[0]*dm4[6 + k1] + 2*c[1]*dm4[15 + k1] + 2*c[2]*dm4[24 + k1] + c[4]*dm4[42 + k1] + 2*c[5]*dm4[51 + k1] + c[8]*dm4[78 + k1] + d[0]*dm5[6 + k1] + 3*d[1]*dm5[15 + k1] + 3*d[2]*dm5[24 + k1] + 3*d[4]*dm5[42 + k1] + 6*d[5]*dm5[51 + k1] + 3*d[8]*dm5[78 + k1] + d[13]*dm5[123 + k1] + 3*d[14]*dm5[132 + k1] + 3*d[17]*dm5[159 + k1] + d[26]*dm5[240 + k1]);
    rhs(7) = -(a*dm2[12 + k1] + b[0]*dm3[12 + k1] + b[1]*dm3[39 + k1] + b[2]*dm3[42 + k1] + c[0]*dm4[12 + k1] + 2*c[1]*dm4[39 + k1] + 2*c[2]*dm4[42 + k1] + c[4]*dm4[120 + k1] + 2*c[5]*dm4[123 + k1] + c[8]*dm4[132 + k1] + d[0]*dm5[12 + k1] + 3*d[1]*dm5[39 + k1] + 3*d[2]*dm5[42 + k1] + 3*d[4]*dm5[120 + k1] + 6*d[5]*dm5[123 + k1] + 3*d[8]*dm5[132 + k1] + d[13]*dm5[363 + k1] + 3*d[14]*dm5[366 + k1] + 3*d[17]*dm5[375 + k1] + d[26]*dm5[402 + k1]);
    rhs(8) = -(a*dm2[15 + k1] + b[0]*dm3[15 + k1] + b[1]*dm3[42 + k1] + b[2]*dm3[51 + k1] + c[0]*dm4[15 + k1] + 2*c[1]*dm4[42 + k1] + 2*c[2]*dm4[51 + k1] + c[4]*dm4[123 + k1] + 2*c[5]*dm4[132 + k1] + c[8]*dm4[159 + k1] + d[0]*dm5[15 + k1] + 3*d[1]*dm5[42 + k1] + 3*d[2]*dm5[51 + k1] + 3*d[4]*dm5[123 + k1] + 6*d[5]*dm5[132 + k1] + 3*d[8]*dm5[159 + k1] + d[13]*dm5[366 + k1] + 3*d[14]*dm5[375 + k1] + 3*d[17]*dm5[402 + k1] + d[26]*dm5[483 + k1]);
    rhs(9) = -(a*dm2[24 + k1] + b[0]*dm3[24 + k1] + b[1]*dm3[51 + k1] + b[2]*dm3[78 + k1] + c[0]*dm4[24 + k1] + 2*c[1]*dm4[51 + k1] + 2*c[2]*dm4[78 + k1] + c[4]*dm4[132 + k1] + 2*c[5]*dm4[159 + k1] + c[8]*dm4[240 + k1] + d[0]*dm5[24 + k1] + 3*d[1]*dm5[51 + k1] + 3*d[2]*dm5[78 + k1] + 3*d[4]*dm5[132 + k1] + 6*d[5]*dm5[159 + k1] + 3*d[8]*dm5[240 + k1] + d[13]*dm5[375 + k1] + 3*d[14]*dm5[402 + k1] + 3*d[17]*dm5[483 + k1] + d[26]*dm5[726 + k1]);
    rhs(10) = -(a*dm3[k1] + b[0]*dm4[k1] + b[1]*dm4[3 + k1] + b[2]*dm4[6 + k1] + c[0]*dm5[k1] + 2*c[1]*dm5[3 + k1] + 2*c[2]*dm5[6 + k1] + c[4]*dm5[12 + k1] + 2*c[5]*dm5[15 + k1] + c[8]*dm5[24 + k1] + d[0]*dm6[k1] + 3*d[1]*dm6[3 + k1] + 3*d[2]*dm6[6 + k1] + 3*d[4]*dm6[12 + k1] + 6*d[5]*dm6[15 + k1] + 3*d[8]*dm6[24 + k1] + d[13]*dm6[39 + k1] + 3*d[14]*dm6[42 + k1] + 3*d[17]*dm6[51 + k1] + d[26]*dm6[78 + k1]);
    rhs(11) = -(a*dm3[3 + k1] + b[0]*dm4[3 + k1] + b[1]*dm4[12 + k1] + b[2]*dm4[15 + k1] + c[0]*dm5[3 + k1] + 2*c[1]*dm5[12 + k1] + 2*c[2]*dm5[15 + k1] + c[4]*dm5[39 + k1] + 2*c[5]*dm5[42 + k1] + c[8]*dm5[51 + k1] + d[0]*dm6[3 + k1] + 3*d[1]*dm6[12 + k1] + 3*d[2]*dm6[15 + k1] + 3*d[4]*dm6[39 + k1] + 6*d[5]*dm6[42 + k1] + 3*d[8]*dm6[51 + k1] + d[13]*dm6[120 + k1] + 3*d[14]*dm6[123 + k1] + 3*d[17]*dm6[132 + k1] + d[26]*dm6[159 + k1]);
    rhs(12) = -(a*dm3[6 + k1] + b[0]*dm4[6 + k1] + b[1]*dm4[15 + k1] + b[2]*dm4[24 + k1] + c[0]*dm5[6 + k1] + 2*c[1]*dm5[15 + k1] + 2*c[2]*dm5[24 + k1] + c[4]*dm5[42 + k1] + 2*c[5]*dm5[51 + k1] + c[8]*dm5[78 + k1] + d[0]*dm6[6 + k1] + 3*d[1]*dm6[15 + k1] + 3*d[2]*dm6[24 + k1] + 3*d[4]*dm6[42 + k1] + 6*d[5]*dm6[51 + k1] + 3*d[8]*dm6[78 + k1] + d[13]*dm6[123 + k1] + 3*d[14]*dm6[132 + k1] + 3*d[17]*dm6[159 + k1] + d[26]*dm6[240 + k1]);
    rhs(13) = -(a*dm3[12 + k1] + b[0]*dm4[12 + k1] + b[1]*dm4[39 + k1] + b[2]*dm4[42 + k1] + c[0]*dm5[12 + k1] + 2*c[1]*dm5[39 + k1] + 2*c[2]*dm5[42 + k1] + c[4]*dm5[120 + k1] + 2*c[5]*dm5[123 + k1] + c[8]*dm5[132 + k1] + d[0]*dm6[12 + k1] + 3*d[1]*dm6[39 + k1] + 3*d[2]*dm6[42 + k1] + 3*d[4]*dm6[120 + k1] + 6*d[5]*dm6[123 + k1] + 3*d[8]*dm6[132 + k1] + d[13]*dm6[363 + k1] + 3*d[14]*dm6[366 + k1] + 3*d[17]*dm6[375 + k1] + d[26]*dm6[402 + k1]);
    rhs(14) = -(a*dm3[15 + k1] + b[0]*dm4[15 + k1] + b[1]*dm4[42 + k1] + b[2]*dm4[51 + k1] + c[0]*dm5[15 + k1] + 2*c[1]*dm5[42 + k1] + 2*c[2]*dm5[51 + k1] + c[4]*dm5[123 + k1] + 2*c[5]*dm5[132 + k1] + c[8]*dm5[159 + k1] + d[0]*dm6[15 + k1] + 3*d[1]*dm6[42 + k1] + 3*d[2]*dm6[51 + k1] + 3*d[4]*dm6[123 + k1] + 6*d[5]*dm6[132 + k1] + 3*d[8]*dm6[159 + k1] + d[13]*dm6[366 + k1] + 3*d[14]*dm6[375 + k1] + 3*d[17]*dm6[402 + k1] + d[26]*dm6[483 + k1]);
    rhs(15) = -(a*dm3[24 + k1] + b[0]*dm4[24 + k1] + b[1]*dm4[51 + k1] + b[2]*dm4[78 + k1] + c[0]*dm5[24 + k1] + 2*c[1]*dm5[51 + k1] + 2*c[2]*dm5[78 + k1] + c[4]*dm5[132 + k1] + 2*c[5]*dm5[159 + k1] + c[8]*dm5[240 + k1] + d[0]*dm6[24 + k1] + 3*d[1]*dm6[51 + k1] + 3*d[2]*dm6[78 + k1] + 3*d[4]*dm6[132 + k1] + 6*d[5]*dm6[159 + k1] + 3*d[8]*dm6[240 + k1] + d[13]*dm6[375 + k1] + 3*d[14]*dm6[402 + k1] + 3*d[17]*dm6[483 + k1] + d[26]*dm6[726 + k1]);
    rhs(16) = -(a*dm3[39 + k1] + b[0]*dm4[39 + k1] + b[1]*dm4[120 + k1] + b[2]*dm4[123 + k1] + c[0]*dm5[39 + k1] + 2*c[1]*dm5[120 + k1] + 2*c[2]*dm5[123 + k1] + c[4]*dm5[363 + k1] + 2*c[5]*dm5[366 + k1] + c[8]*dm5[375 + k1] + d[0]*dm6[39 + k1] + 3*d[1]*dm6[120 + k1] + 3*d[2]*dm6[123 + k1] + 3*d[4]*dm6[363 + k1] + 6*d[5]*dm6[366 + k1] + 3*d[8]*dm6[375 + k1] + d[13]*dm6[1092 + k1] + 3*d[14]*dm6[1095 + k1] + 3*d[17]*dm6[1104 + k1] + d[26]*dm6[1131 + k1]);
    rhs(17) = -(a*dm3[42 + k1] + b[0]*dm4[42 + k1] + b[1]*dm4[123 + k1] + b[2]*dm4[132 + k1] + c[0]*dm5[42 + k1] + 2*c[1]*dm5[123 + k1] + 2*c[2]*dm5[132 + k1] + c[4]*dm5[366 + k1] + 2*c[5]*dm5[375 + k1] + c[8]*dm5[402 + k1] + d[0]*dm6[42 + k1] + 3*d[1]*dm6[123 + k1] + 3*d[2]*dm6[132 + k1] + 3*d[4]*dm6[366 + k1] + 6*d[5]*dm6[375 + k1] + 3*d[8]*dm6[402 + k1] + d[13]*dm6[1095 + k1] + 3*d[14]*dm6[1104 + k1] + 3*d[17]*dm6[1131 + k1] + d[26]*dm6[1212 + k1]);
    rhs(18) = -(a*dm3[51 + k1] + b[0]*dm4[51 + k1] + b[1]*dm4[132 + k1] + b[2]*dm4[159 + k1] + c[0]*dm5[51 + k1] + 2*c[1]*dm5[132 + k1] + 2*c[2]*dm5[159 + k1] + c[4]*dm5[375 + k1] + 2*c[5]*dm5[402 + k1] + c[8]*dm5[483 + k1] + d[0]*dm6[51 + k1] + 3*d[1]*dm6[132 + k1] + 3*d[2]*dm6[159 + k1] + 3*d[4]*dm6[375 + k1] + 6*d[5]*dm6[402 + k1] + 3*d[8]*dm6[483 + k1] + d[13]*dm6[1104 + k1] + 3*d[14]*dm6[1131 + k1] + 3*d[17]*dm6[1212 + k1] + d[26]*dm6[1455 + k1]);
    rhs(19) = -(a*dm3[78 + k1] + b[0]*dm4[78 + k1] + b[1]*dm4[159 + k1] + b[2]*dm4[240 + k1] + c[0]*dm5[78 + k1] + 2*c[1]*dm5[159 + k1] + 2*c[2]*dm5[240 + k1] + c[4]*dm5[402 + k1] + 2*c[5]*dm5[483 + k1] + c[8]*dm5[726 + k1] + d[0]*dm6[78 + k1] + 3*d[1]*dm6[159 + k1] + 3*d[2]*dm6[240 + k1] + 3*d[4]*dm6[402 + k1] + 6*d[5]*dm6[483 + k1] + 3*d[8]*dm6[726 + k1] + d[13]*dm6[1131 + k1] + 3*d[14]*dm6[1212 + k1] + 3*d[17]*dm6[1455 + k1] + d[26]*dm6[2184 + k1]);
 
    lhs = solver.solve(rhs);

    da[k1] = lhs(0);
    db[k1] = lhs(1);
    db[3 + k1] = lhs(2);
    db[6 + k1] = lhs(3);
    dc[k1] = lhs(4);
    dc[3 + k1] = lhs(5);
    dc[6 + k1] = lhs(6);
    dc[12 + k1] = lhs(7);
    dc[15 + k1] = lhs(8);
    dc[24 + k1] = lhs(9);
    dd[k1] = lhs(10);
    dd[3 + k1] = lhs(11);
    dd[6 + k1] = lhs(12);
    dd[12 + k1] = lhs(13);
    dd[15 + k1] = lhs(14);
    dd[24 + k1] = lhs(15);
    dd[39 + k1] = lhs(16);
    dd[42 + k1] = lhs(17);
    dd[51 + k1] = lhs(18);
    dd[78 + k1] = lhs(19);
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[3*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[3*k1 + k2] + b[1]*ddm1[9 + 3*k1 + k2] + b[2]*ddm1[18 + 3*k1 + k2] + c[0]*ddm2[3*k1 + k2] + 2*c[1]*ddm2[9 + 3*k1 + k2] + 2*c[2]*ddm2[18 + 3*k1 + k2] + c[4]*ddm2[36 + 3*k1 + k2] + 2*c[5]*ddm2[45 + 3*k1 + k2] + c[8]*ddm2[72 + 3*k1 + k2] + d[0]*ddm3[3*k1 + k2] + 3*d[1]*ddm3[9 + 3*k1 + k2] + 3*d[2]*ddm3[18 + 3*k1 + k2] + 3*d[4]*ddm3[36 + 3*k1 + k2] + 6*d[5]*ddm3[45 + 3*k1 + k2] + 3*d[8]*ddm3[72 + 3*k1 + k2] + d[13]*ddm3[117 + 3*k1 + k2] + 3*d[14]*ddm3[126 + 3*k1 + k2] + 3*d[17]*ddm3[153 + 3*k1 + k2] + d[26]*ddm3[234 + 3*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[3 + k2]*dm1[3 + k1] + db[3 + k1]*dm1[3 + k2] + db[6 + k2]*dm1[6 + k1] + db[6 + k1]*dm1[6 + k2] + dc[k2]*dm2[k1] + dc[k1]*dm2[k2] + 2*dc[3 + k2]*dm2[3 + k1] + 2*dc[3 + k1]*dm2[3 + k2] + 2*dc[6 + k2]*dm2[6 + k1] + 2*dc[6 + k1]*dm2[6 + k2] + dc[12 + k2]*dm2[12 + k1] + dc[12 + k1]*dm2[12 + k2] + 2*dc[15 + k2]*dm2[15 + k1] + 2*dc[15 + k1]*dm2[15 + k2] + dc[24 + k2]*dm2[24 + k1] + dc[24 + k1]*dm2[24 + k2] + dd[k2]*dm3[k1] + dd[k1]*dm3[k2] + 3*dd[3 + k2]*dm3[3 + k1] + 3*dd[3 + k1]*dm3[3 + k2] + 3*dd[6 + k2]*dm3[6 + k1] + 3*dd[6 + k1]*dm3[6 + k2] + 3*dd[12 + k2]*dm3[12 + k1] + 3*dd[12 + k1]*dm3[12 + k2] + 6*dd[15 + k2]*dm3[15 + k1] + 6*dd[15 + k1]*dm3[15 + k2] + 3*dd[24 + k2]*dm3[24 + k1] + 3*dd[24 + k1]*dm3[24 + k2] + dd[39 + k2]*dm3[39 + k1] + dd[39 + k1]*dm3[39 + k2] + 3*dd[42 + k2]*dm3[42 + k1] + 3*dd[42 + k1]*dm3[42 + k2] + 3*dd[51 + k2]*dm3[51 + k1] + 3*dd[51 + k1]*dm3[51 + k2] + dd[78 + k2]*dm3[78 + k1] + dd[78 + k1]*dm3[78 + k2]);
      rhs(1) = -(a*ddm1[3*k1 + k2] + b[0]*ddm2[3*k1 + k2] + b[1]*ddm2[9 + 3*k1 + k2] + b[2]*ddm2[18 + 3*k1 + k2] + c[0]*ddm3[3*k1 + k2] + 2*c[1]*ddm3[9 + 3*k1 + k2] + 2*c[2]*ddm3[18 + 3*k1 + k2] + c[4]*ddm3[36 + 3*k1 + k2] + 2*c[5]*ddm3[45 + 3*k1 + k2] + c[8]*ddm3[72 + 3*k1 + k2] + d[0]*ddm4[3*k1 + k2] + 3*d[1]*ddm4[9 + 3*k1 + k2] + 3*d[2]*ddm4[18 + 3*k1 + k2] + 3*d[4]*ddm4[36 + 3*k1 + k2] + 6*d[5]*ddm4[45 + 3*k1 + k2] + 3*d[8]*ddm4[72 + 3*k1 + k2] + d[13]*ddm4[117 + 3*k1 + k2] + 3*d[14]*ddm4[126 + 3*k1 + k2] + 3*d[17]*ddm4[153 + 3*k1 + k2] + d[26]*ddm4[234 + 3*k1 + k2] + da[k2]*dm1[k1] + da[k1]*dm1[k2] + db[k2]*dm2[k1] + db[k1]*dm2[k2] + db[3 + k2]*dm2[3 + k1] + db[3 + k1]*dm2[3 + k2] + db[6 + k2]*dm2[6 + k1] + db[6 + k1]*dm2[6 + k2] + dc[k2]*dm3[k1] + dc[k1]*dm3[k2] + 2*dc[3 + k2]*dm3[3 + k1] + 2*dc[3 + k1]*dm3[3 + k2] + 2*dc[6 + k2]*dm3[6 + k1] + 2*dc[6 + k1]*dm3[6 + k2] + dc[12 + k2]*dm3[12 + k1] + dc[12 + k1]*dm3[12 + k2] + 2*dc[15 + k2]*dm3[15 + k1] + 2*dc[15 + k1]*dm3[15 + k2] + dc[24 + k2]*dm3[24 + k1] + dc[24 + k1]*dm3[24 + k2] + dd[k2]*dm4[k1] + dd[k1]*dm4[k2] + 3*dd[3 + k2]*dm4[3 + k1] + 3*dd[3 + k1]*dm4[3 + k2] + 3*dd[6 + k2]*dm4[6 + k1] + 3*dd[6 + k1]*dm4[6 + k2] + 3*dd[12 + k2]*dm4[12 + k1] + 3*dd[12 + k1]*dm4[12 + k2] + 6*dd[15 + k2]*dm4[15 + k1] + 6*dd[15 + k1]*dm4[15 + k2] + 3*dd[24 + k2]*dm4[24 + k1] + 3*dd[24 + k1]*dm4[24 + k2] + dd[39 + k2]*dm4[39 + k1] + dd[39 + k1]*dm4[39 + k2] + 3*dd[42 + k2]*dm4[42 + k1] + 3*dd[42 + k1]*dm4[42 + k2] + 3*dd[51 + k2]*dm4[51 + k1] + 3*dd[51 + k1]*dm4[51 + k2] + dd[78 + k2]*dm4[78 + k1] + dd[78 + k1]*dm4[78 + k2]);
      rhs(2) = -(a*ddm1[9 + 3*k1 + k2] + b[0]*ddm2[9 + 3*k1 + k2] + b[1]*ddm2[36 + 3*k1 + k2] + b[2]*ddm2[45 + 3*k1 + k2] + c[0]*ddm3[9 + 3*k1 + k2] + 2*c[1]*ddm3[36 + 3*k1 + k2] + 2*c[2]*ddm3[45 + 3*k1 + k2] + c[4]*ddm3[117 + 3*k1 + k2] + 2*c[5]*ddm3[126 + 3*k1 + k2] + c[8]*ddm3[153 + 3*k1 + k2] + d[0]*ddm4[9 + 3*k1 + k2] + 3*d[1]*ddm4[36 + 3*k1 + k2] + 3*d[2]*ddm4[45 + 3*k1 + k2] + 3*d[4]*ddm4[117 + 3*k1 + k2] + 6*d[5]*ddm4[126 + 3*k1 + k2] + 3*d[8]*ddm4[153 + 3*k1 + k2] + d[13]*ddm4[360 + 3*k1 + k2] + 3*d[14]*ddm4[369 + 3*k1 + k2] + 3*d[17]*ddm4[396 + 3*k1 + k2] + d[26]*ddm4[477 + 3*k1 + k2] + da[k2]*dm1[3 + k1] + da[k1]*dm1[3 + k2] + db[k2]*dm2[3 + k1] + db[k1]*dm2[3 + k2] + db[3 + k2]*dm2[12 + k1] + db[3 + k1]*dm2[12 + k2] + db[6 + k2]*dm2[15 + k1] + db[6 + k1]*dm2[15 + k2] + dc[k2]*dm3[3 + k1] + dc[k1]*dm3[3 + k2] + 2*dc[3 + k2]*dm3[12 + k1] + 2*dc[3 + k1]*dm3[12 + k2] + 2*dc[6 + k2]*dm3[15 + k1] + 2*dc[6 + k1]*dm3[15 + k2] + dc[12 + k2]*dm3[39 + k1] + dc[12 + k1]*dm3[39 + k2] + 2*dc[15 + k2]*dm3[42 + k1] + 2*dc[15 + k1]*dm3[42 + k2] + dc[24 + k2]*dm3[51 + k1] + dc[24 + k1]*dm3[51 + k2] + dd[k2]*dm4[3 + k1] + dd[k1]*dm4[3 + k2] + 3*dd[3 + k2]*dm4[12 + k1] + 3*dd[3 + k1]*dm4[12 + k2] + 3*dd[6 + k2]*dm4[15 + k1] + 3*dd[6 + k1]*dm4[15 + k2] + 3*dd[12 + k2]*dm4[39 + k1] + 3*dd[12 + k1]*dm4[39 + k2] + 6*dd[15 + k2]*dm4[42 + k1] + 6*dd[15 + k1]*dm4[42 + k2] + 3*dd[24 + k2]*dm4[51 + k1] + 3*dd[24 + k1]*dm4[51 + k2] + dd[39 + k2]*dm4[120 + k1] + dd[39 + k1]*dm4[120 + k2] + 3*dd[42 + k2]*dm4[123 + k1] + 3*dd[42 + k1]*dm4[123 + k2] + 3*dd[51 + k2]*dm4[132 + k1] + 3*dd[51 + k1]*dm4[132 + k2] + dd[78 + k2]*dm4[159 + k1] + dd[78 + k1]*dm4[159 + k2]);
      rhs(3) = -(a*ddm1[18 + 3*k1 + k2] + b[0]*ddm2[18 + 3*k1 + k2] + b[1]*ddm2[45 + 3*k1 + k2] + b[2]*ddm2[72 + 3*k1 + k2] + c[0]*ddm3[18 + 3*k1 + k2] + 2*c[1]*ddm3[45 + 3*k1 + k2] + 2*c[2]*ddm3[72 + 3*k1 + k2] + c[4]*ddm3[126 + 3*k1 + k2] + 2*c[5]*ddm3[153 + 3*k1 + k2] + c[8]*ddm3[234 + 3*k1 + k2] + d[0]*ddm4[18 + 3*k1 + k2] + 3*d[1]*ddm4[45 + 3*k1 + k2] + 3*d[2]*ddm4[72 + 3*k1 + k2] + 3*d[4]*ddm4[126 + 3*k1 + k2] + 6*d[5]*ddm4[153 + 3*k1 + k2] + 3*d[8]*ddm4[234 + 3*k1 + k2] + d[13]*ddm4[369 + 3*k1 + k2] + 3*d[14]*ddm4[396 + 3*k1 + k2] + 3*d[17]*ddm4[477 + 3*k1 + k2] + d[26]*ddm4[720 + 3*k1 + k2] + da[k2]*dm1[6 + k1] + da[k1]*dm1[6 + k2] + db[k2]*dm2[6 + k1] + db[k1]*dm2[6 + k2] + db[3 + k2]*dm2[15 + k1] + db[3 + k1]*dm2[15 + k2] + db[6 + k2]*dm2[24 + k1] + db[6 + k1]*dm2[24 + k2] + dc[k2]*dm3[6 + k1] + dc[k1]*dm3[6 + k2] + 2*dc[3 + k2]*dm3[15 + k1] + 2*dc[3 + k1]*dm3[15 + k2] + 2*dc[6 + k2]*dm3[24 + k1] + 2*dc[6 + k1]*dm3[24 + k2] + dc[12 + k2]*dm3[42 + k1] + dc[12 + k1]*dm3[42 + k2] + 2*dc[15 + k2]*dm3[51 + k1] + 2*dc[15 + k1]*dm3[51 + k2] + dc[24 + k2]*dm3[78 + k1] + dc[24 + k1]*dm3[78 + k2] + dd[k2]*dm4[6 + k1] + dd[k1]*dm4[6 + k2] + 3*dd[3 + k2]*dm4[15 + k1] + 3*dd[3 + k1]*dm4[15 + k2] + 3*dd[6 + k2]*dm4[24 + k1] + 3*dd[6 + k1]*dm4[24 + k2] + 3*dd[12 + k2]*dm4[42 + k1] + 3*dd[12 + k1]*dm4[42 + k2] + 6*dd[15 + k2]*dm4[51 + k1] + 6*dd[15 + k1]*dm4[51 + k2] + 3*dd[24 + k2]*dm4[78 + k1] + 3*dd[24 + k1]*dm4[78 + k2] + dd[39 + k2]*dm4[123 + k1] + dd[39 + k1]*dm4[123 + k2] + 3*dd[42 + k2]*dm4[132 + k1] + 3*dd[42 + k1]*dm4[132 + k2] + 3*dd[51 + k2]*dm4[159 + k1] + 3*dd[51 + k1]*dm4[159 + k2] + dd[78 + k2]*dm4[240 + k1] + dd[78 + k1]*dm4[240 + k2]);
      rhs(4) = -(a*ddm2[3*k1 + k2] + b[0]*ddm3[3*k1 + k2] + b[1]*ddm3[9 + 3*k1 + k2] + b[2]*ddm3[18 + 3*k1 + k2] + c[0]*ddm4[3*k1 + k2] + 2*c[1]*ddm4[9 + 3*k1 + k2] + 2*c[2]*ddm4[18 + 3*k1 + k2] + c[4]*ddm4[36 + 3*k1 + k2] + 2*c[5]*ddm4[45 + 3*k1 + k2] + c[8]*ddm4[72 + 3*k1 + k2] + d[0]*ddm5[3*k1 + k2] + 3*d[1]*ddm5[9 + 3*k1 + k2] + 3*d[2]*ddm5[18 + 3*k1 + k2] + 3*d[4]*ddm5[36 + 3*k1 + k2] + 6*d[5]*ddm5[45 + 3*k1 + k2] + 3*d[8]*ddm5[72 + 3*k1 + k2] + d[13]*ddm5[117 + 3*k1 + k2] + 3*d[14]*ddm5[126 + 3*k1 + k2] + 3*d[17]*ddm5[153 + 3*k1 + k2] + d[26]*ddm5[234 + 3*k1 + k2] + da[k2]*dm2[k1] + da[k1]*dm2[k2] + db[k2]*dm3[k1] + db[k1]*dm3[k2] + db[3 + k2]*dm3[3 + k1] + db[3 + k1]*dm3[3 + k2] + db[6 + k2]*dm3[6 + k1] + db[6 + k1]*dm3[6 + k2] + dc[k2]*dm4[k1] + dc[k1]*dm4[k2] + 2*dc[3 + k2]*dm4[3 + k1] + 2*dc[3 + k1]*dm4[3 + k2] + 2*dc[6 + k2]*dm4[6 + k1] + 2*dc[6 + k1]*dm4[6 + k2] + dc[12 + k2]*dm4[12 + k1] + dc[12 + k1]*dm4[12 + k2] + 2*dc[15 + k2]*dm4[15 + k1] + 2*dc[15 + k1]*dm4[15 + k2] + dc[24 + k2]*dm4[24 + k1] + dc[24 + k1]*dm4[24 + k2] + dd[k2]*dm5[k1] + dd[k1]*dm5[k2] + 3*dd[3 + k2]*dm5[3 + k1] + 3*dd[3 + k1]*dm5[3 + k2] + 3*dd[6 + k2]*dm5[6 + k1] + 3*dd[6 + k1]*dm5[6 + k2] + 3*dd[12 + k2]*dm5[12 + k1] + 3*dd[12 + k1]*dm5[12 + k2] + 6*dd[15 + k2]*dm5[15 + k1] + 6*dd[15 + k1]*dm5[15 + k2] + 3*dd[24 + k2]*dm5[24 + k1] + 3*dd[24 + k1]*dm5[24 + k2] + dd[39 + k2]*dm5[39 + k1] + dd[39 + k1]*dm5[39 + k2] + 3*dd[42 + k2]*dm5[42 + k1] + 3*dd[42 + k1]*dm5[42 + k2] + 3*dd[51 + k2]*dm5[51 + k1] + 3*dd[51 + k1]*dm5[51 + k2] + dd[78 + k2]*dm5[78 + k1] + dd[78 + k1]*dm5[78 + k2]);
      rhs(5) = -(a*ddm2[9 + 3*k1 + k2] + b[0]*ddm3[9 + 3*k1 + k2] + b[1]*ddm3[36 + 3*k1 + k2] + b[2]*ddm3[45 + 3*k1 + k2] + c[0]*ddm4[9 + 3*k1 + k2] + 2*c[1]*ddm4[36 + 3*k1 + k2] + 2*c[2]*ddm4[45 + 3*k1 + k2] + c[4]*ddm4[117 + 3*k1 + k2] + 2*c[5]*ddm4[126 + 3*k1 + k2] + c[8]*ddm4[153 + 3*k1 + k2] + d[0]*ddm5[9 + 3*k1 + k2] + 3*d[1]*ddm5[36 + 3*k1 + k2] + 3*d[2]*ddm5[45 + 3*k1 + k2] + 3*d[4]*ddm5[117 + 3*k1 + k2] + 6*d[5]*ddm5[126 + 3*k1 + k2] + 3*d[8]*ddm5[153 + 3*k1 + k2] + d[13]*ddm5[360 + 3*k1 + k2] + 3*d[14]*ddm5[369 + 3*k1 + k2] + 3*d[17]*ddm5[396 + 3*k1 + k2] + d[26]*ddm5[477 + 3*k1 + k2] + da[k2]*dm2[3 + k1] + da[k1]*dm2[3 + k2] + db[k2]*dm3[3 + k1] + db[k1]*dm3[3 + k2] + db[3 + k2]*dm3[12 + k1] + db[3 + k1]*dm3[12 + k2] + db[6 + k2]*dm3[15 + k1] + db[6 + k1]*dm3[15 + k2] + dc[k2]*dm4[3 + k1] + dc[k1]*dm4[3 + k2] + 2*dc[3 + k2]*dm4[12 + k1] + 2*dc[3 + k1]*dm4[12 + k2] + 2*dc[6 + k2]*dm4[15 + k1] + 2*dc[6 + k1]*dm4[15 + k2] + dc[12 + k2]*dm4[39 + k1] + dc[12 + k1]*dm4[39 + k2] + 2*dc[15 + k2]*dm4[42 + k1] + 2*dc[15 + k1]*dm4[42 + k2] + dc[24 + k2]*dm4[51 + k1] + dc[24 + k1]*dm4[51 + k2] + dd[k2]*dm5[3 + k1] + dd[k1]*dm5[3 + k2] + 3*dd[3 + k2]*dm5[12 + k1] + 3*dd[3 + k1]*dm5[12 + k2] + 3*dd[6 + k2]*dm5[15 + k1] + 3*dd[6 + k1]*dm5[15 + k2] + 3*dd[12 + k2]*dm5[39 + k1] + 3*dd[12 + k1]*dm5[39 + k2] + 6*dd[15 + k2]*dm5[42 + k1] + 6*dd[15 + k1]*dm5[42 + k2] + 3*dd[24 + k2]*dm5[51 + k1] + 3*dd[24 + k1]*dm5[51 + k2] + dd[39 + k2]*dm5[120 + k1] + dd[39 + k1]*dm5[120 + k2] + 3*dd[42 + k2]*dm5[123 + k1] + 3*dd[42 + k1]*dm5[123 + k2] + 3*dd[51 + k2]*dm5[132 + k1] + 3*dd[51 + k1]*dm5[132 + k2] + dd[78 + k2]*dm5[159 + k1] + dd[78 + k1]*dm5[159 + k2]);
      rhs(6) = -(a*ddm2[18 + 3*k1 + k2] + b[0]*ddm3[18 + 3*k1 + k2] + b[1]*ddm3[45 + 3*k1 + k2] + b[2]*ddm3[72 + 3*k1 + k2] + c[0]*ddm4[18 + 3*k1 + k2] + 2*c[1]*ddm4[45 + 3*k1 + k2] + 2*c[2]*ddm4[72 + 3*k1 + k2] + c[4]*ddm4[126 + 3*k1 + k2] + 2*c[5]*ddm4[153 + 3*k1 + k2] + c[8]*ddm4[234 + 3*k1 + k2] + d[0]*ddm5[18 + 3*k1 + k2] + 3*d[1]*ddm5[45 + 3*k1 + k2] + 3*d[2]*ddm5[72 + 3*k1 + k2] + 3*d[4]*ddm5[126 + 3*k1 + k2] + 6*d[5]*ddm5[153 + 3*k1 + k2] + 3*d[8]*ddm5[234 + 3*k1 + k2] + d[13]*ddm5[369 + 3*k1 + k2] + 3*d[14]*ddm5[396 + 3*k1 + k2] + 3*d[17]*ddm5[477 + 3*k1 + k2] + d[26]*ddm5[720 + 3*k1 + k2] + da[k2]*dm2[6 + k1] + da[k1]*dm2[6 + k2] + db[k2]*dm3[6 + k1] + db[k1]*dm3[6 + k2] + db[3 + k2]*dm3[15 + k1] + db[3 + k1]*dm3[15 + k2] + db[6 + k2]*dm3[24 + k1] + db[6 + k1]*dm3[24 + k2] + dc[k2]*dm4[6 + k1] + dc[k1]*dm4[6 + k2] + 2*dc[3 + k2]*dm4[15 + k1] + 2*dc[3 + k1]*dm4[15 + k2] + 2*dc[6 + k2]*dm4[24 + k1] + 2*dc[6 + k1]*dm4[24 + k2] + dc[12 + k2]*dm4[42 + k1] + dc[12 + k1]*dm4[42 + k2] + 2*dc[15 + k2]*dm4[51 + k1] + 2*dc[15 + k1]*dm4[51 + k2] + dc[24 + k2]*dm4[78 + k1] + dc[24 + k1]*dm4[78 + k2] + dd[k2]*dm5[6 + k1] + dd[k1]*dm5[6 + k2] + 3*dd[3 + k2]*dm5[15 + k1] + 3*dd[3 + k1]*dm5[15 + k2] + 3*dd[6 + k2]*dm5[24 + k1] + 3*dd[6 + k1]*dm5[24 + k2] + 3*dd[12 + k2]*dm5[42 + k1] + 3*dd[12 + k1]*dm5[42 + k2] + 6*dd[15 + k2]*dm5[51 + k1] + 6*dd[15 + k1]*dm5[51 + k2] + 3*dd[24 + k2]*dm5[78 + k1] + 3*dd[24 + k1]*dm5[78 + k2] + dd[39 + k2]*dm5[123 + k1] + dd[39 + k1]*dm5[123 + k2] + 3*dd[42 + k2]*dm5[132 + k1] + 3*dd[42 + k1]*dm5[132 + k2] + 3*dd[51 + k2]*dm5[159 + k1] + 3*dd[51 + k1]*dm5[159 + k2] + dd[78 + k2]*dm5[240 + k1] + dd[78 + k1]*dm5[240 + k2]);
      rhs(7) = -(a*ddm2[36 + 3*k1 + k2] + b[0]*ddm3[36 + 3*k1 + k2] + b[1]*ddm3[117 + 3*k1 + k2] + b[2]*ddm3[126 + 3*k1 + k2] + c[0]*ddm4[36 + 3*k1 + k2] + 2*c[1]*ddm4[117 + 3*k1 + k2] + 2*c[2]*ddm4[126 + 3*k1 + k2] + c[4]*ddm4[360 + 3*k1 + k2] + 2*c[5]*ddm4[369 + 3*k1 + k2] + c[8]*ddm4[396 + 3*k1 + k2] + d[0]*ddm5[36 + 3*k1 + k2] + 3*d[1]*ddm5[117 + 3*k1 + k2] + 3*d[2]*ddm5[126 + 3*k1 + k2] + 3*d[4]*ddm5[360 + 3*k1 + k2] + 6*d[5]*ddm5[369 + 3*k1 + k2] + 3*d[8]*ddm5[396 + 3*k1 + k2] + d[13]*ddm5[1089 + 3*k1 + k2] + 3*d[14]*ddm5[1098 + 3*k1 + k2] + 3*d[17]*ddm5[1125 + 3*k1 + k2] + d[26]*ddm5[1206 + 3*k1 + k2] + da[k2]*dm2[12 + k1] + da[k1]*dm2[12 + k2] + db[k2]*dm3[12 + k1] + db[k1]*dm3[12 + k2] + db[3 + k2]*dm3[39 + k1] + db[3 + k1]*dm3[39 + k2] + db[6 + k2]*dm3[42 + k1] + db[6 + k1]*dm3[42 + k2] + dc[k2]*dm4[12 + k1] + dc[k1]*dm4[12 + k2] + 2*dc[3 + k2]*dm4[39 + k1] + 2*dc[3 + k1]*dm4[39 + k2] + 2*dc[6 + k2]*dm4[42 + k1] + 2*dc[6 + k1]*dm4[42 + k2] + dc[12 + k2]*dm4[120 + k1] + dc[12 + k1]*dm4[120 + k2] + 2*dc[15 + k2]*dm4[123 + k1] + 2*dc[15 + k1]*dm4[123 + k2] + dc[24 + k2]*dm4[132 + k1] + dc[24 + k1]*dm4[132 + k2] + dd[k2]*dm5[12 + k1] + dd[k1]*dm5[12 + k2] + 3*dd[3 + k2]*dm5[39 + k1] + 3*dd[3 + k1]*dm5[39 + k2] + 3*dd[6 + k2]*dm5[42 + k1] + 3*dd[6 + k1]*dm5[42 + k2] + 3*dd[12 + k2]*dm5[120 + k1] + 3*dd[12 + k1]*dm5[120 + k2] + 6*dd[15 + k2]*dm5[123 + k1] + 6*dd[15 + k1]*dm5[123 + k2] + 3*dd[24 + k2]*dm5[132 + k1] + 3*dd[24 + k1]*dm5[132 + k2] + dd[39 + k2]*dm5[363 + k1] + dd[39 + k1]*dm5[363 + k2] + 3*dd[42 + k2]*dm5[366 + k1] + 3*dd[42 + k1]*dm5[366 + k2] + 3*dd[51 + k2]*dm5[375 + k1] + 3*dd[51 + k1]*dm5[375 + k2] + dd[78 + k2]*dm5[402 + k1] + dd[78 + k1]*dm5[402 + k2]);
      rhs(8) = -(a*ddm2[45 + 3*k1 + k2] + b[0]*ddm3[45 + 3*k1 + k2] + b[1]*ddm3[126 + 3*k1 + k2] + b[2]*ddm3[153 + 3*k1 + k2] + c[0]*ddm4[45 + 3*k1 + k2] + 2*c[1]*ddm4[126 + 3*k1 + k2] + 2*c[2]*ddm4[153 + 3*k1 + k2] + c[4]*ddm4[369 + 3*k1 + k2] + 2*c[5]*ddm4[396 + 3*k1 + k2] + c[8]*ddm4[477 + 3*k1 + k2] + d[0]*ddm5[45 + 3*k1 + k2] + 3*d[1]*ddm5[126 + 3*k1 + k2] + 3*d[2]*ddm5[153 + 3*k1 + k2] + 3*d[4]*ddm5[369 + 3*k1 + k2] + 6*d[5]*ddm5[396 + 3*k1 + k2] + 3*d[8]*ddm5[477 + 3*k1 + k2] + d[13]*ddm5[1098 + 3*k1 + k2] + 3*d[14]*ddm5[1125 + 3*k1 + k2] + 3*d[17]*ddm5[1206 + 3*k1 + k2] + d[26]*ddm5[1449 + 3*k1 + k2] + da[k2]*dm2[15 + k1] + da[k1]*dm2[15 + k2] + db[k2]*dm3[15 + k1] + db[k1]*dm3[15 + k2] + db[3 + k2]*dm3[42 + k1] + db[3 + k1]*dm3[42 + k2] + db[6 + k2]*dm3[51 + k1] + db[6 + k1]*dm3[51 + k2] + dc[k2]*dm4[15 + k1] + dc[k1]*dm4[15 + k2] + 2*dc[3 + k2]*dm4[42 + k1] + 2*dc[3 + k1]*dm4[42 + k2] + 2*dc[6 + k2]*dm4[51 + k1] + 2*dc[6 + k1]*dm4[51 + k2] + dc[12 + k2]*dm4[123 + k1] + dc[12 + k1]*dm4[123 + k2] + 2*dc[15 + k2]*dm4[132 + k1] + 2*dc[15 + k1]*dm4[132 + k2] + dc[24 + k2]*dm4[159 + k1] + dc[24 + k1]*dm4[159 + k2] + dd[k2]*dm5[15 + k1] + dd[k1]*dm5[15 + k2] + 3*dd[3 + k2]*dm5[42 + k1] + 3*dd[3 + k1]*dm5[42 + k2] + 3*dd[6 + k2]*dm5[51 + k1] + 3*dd[6 + k1]*dm5[51 + k2] + 3*dd[12 + k2]*dm5[123 + k1] + 3*dd[12 + k1]*dm5[123 + k2] + 6*dd[15 + k2]*dm5[132 + k1] + 6*dd[15 + k1]*dm5[132 + k2] + 3*dd[24 + k2]*dm5[159 + k1] + 3*dd[24 + k1]*dm5[159 + k2] + dd[39 + k2]*dm5[366 + k1] + dd[39 + k1]*dm5[366 + k2] + 3*dd[42 + k2]*dm5[375 + k1] + 3*dd[42 + k1]*dm5[375 + k2] + 3*dd[51 + k2]*dm5[402 + k1] + 3*dd[51 + k1]*dm5[402 + k2] + dd[78 + k2]*dm5[483 + k1] + dd[78 + k1]*dm5[483 + k2]);
      rhs(9) = -(a*ddm2[72 + 3*k1 + k2] + b[0]*ddm3[72 + 3*k1 + k2] + b[1]*ddm3[153 + 3*k1 + k2] + b[2]*ddm3[234 + 3*k1 + k2] + c[0]*ddm4[72 + 3*k1 + k2] + 2*c[1]*ddm4[153 + 3*k1 + k2] + 2*c[2]*ddm4[234 + 3*k1 + k2] + c[4]*ddm4[396 + 3*k1 + k2] + 2*c[5]*ddm4[477 + 3*k1 + k2] + c[8]*ddm4[720 + 3*k1 + k2] + d[0]*ddm5[72 + 3*k1 + k2] + 3*d[1]*ddm5[153 + 3*k1 + k2] + 3*d[2]*ddm5[234 + 3*k1 + k2] + 3*d[4]*ddm5[396 + 3*k1 + k2] + 6*d[5]*ddm5[477 + 3*k1 + k2] + 3*d[8]*ddm5[720 + 3*k1 + k2] + d[13]*ddm5[1125 + 3*k1 + k2] + 3*d[14]*ddm5[1206 + 3*k1 + k2] + 3*d[17]*ddm5[1449 + 3*k1 + k2] + d[26]*ddm5[2178 + 3*k1 + k2] + da[k2]*dm2[24 + k1] + da[k1]*dm2[24 + k2] + db[k2]*dm3[24 + k1] + db[k1]*dm3[24 + k2] + db[3 + k2]*dm3[51 + k1] + db[3 + k1]*dm3[51 + k2] + db[6 + k2]*dm3[78 + k1] + db[6 + k1]*dm3[78 + k2] + dc[k2]*dm4[24 + k1] + dc[k1]*dm4[24 + k2] + 2*dc[3 + k2]*dm4[51 + k1] + 2*dc[3 + k1]*dm4[51 + k2] + 2*dc[6 + k2]*dm4[78 + k1] + 2*dc[6 + k1]*dm4[78 + k2] + dc[12 + k2]*dm4[132 + k1] + dc[12 + k1]*dm4[132 + k2] + 2*dc[15 + k2]*dm4[159 + k1] + 2*dc[15 + k1]*dm4[159 + k2] + dc[24 + k2]*dm4[240 + k1] + dc[24 + k1]*dm4[240 + k2] + dd[k2]*dm5[24 + k1] + dd[k1]*dm5[24 + k2] + 3*dd[3 + k2]*dm5[51 + k1] + 3*dd[3 + k1]*dm5[51 + k2] + 3*dd[6 + k2]*dm5[78 + k1] + 3*dd[6 + k1]*dm5[78 + k2] + 3*dd[12 + k2]*dm5[132 + k1] + 3*dd[12 + k1]*dm5[132 + k2] + 6*dd[15 + k2]*dm5[159 + k1] + 6*dd[15 + k1]*dm5[159 + k2] + 3*dd[24 + k2]*dm5[240 + k1] + 3*dd[24 + k1]*dm5[240 + k2] + dd[39 + k2]*dm5[375 + k1] + dd[39 + k1]*dm5[375 + k2] + 3*dd[42 + k2]*dm5[402 + k1] + 3*dd[42 + k1]*dm5[402 + k2] + 3*dd[51 + k2]*dm5[483 + k1] + 3*dd[51 + k1]*dm5[483 + k2] + dd[78 + k2]*dm5[726 + k1] + dd[78 + k1]*dm5[726 + k2]);
      rhs(10) = -(a*ddm3[3*k1 + k2] + b[0]*ddm4[3*k1 + k2] + b[1]*ddm4[9 + 3*k1 + k2] + b[2]*ddm4[18 + 3*k1 + k2] + c[0]*ddm5[3*k1 + k2] + 2*c[1]*ddm5[9 + 3*k1 + k2] + 2*c[2]*ddm5[18 + 3*k1 + k2] + c[4]*ddm5[36 + 3*k1 + k2] + 2*c[5]*ddm5[45 + 3*k1 + k2] + c[8]*ddm5[72 + 3*k1 + k2] + d[0]*ddm6[3*k1 + k2] + 3*d[1]*ddm6[9 + 3*k1 + k2] + 3*d[2]*ddm6[18 + 3*k1 + k2] + 3*d[4]*ddm6[36 + 3*k1 + k2] + 6*d[5]*ddm6[45 + 3*k1 + k2] + 3*d[8]*ddm6[72 + 3*k1 + k2] + d[13]*ddm6[117 + 3*k1 + k2] + 3*d[14]*ddm6[126 + 3*k1 + k2] + 3*d[17]*ddm6[153 + 3*k1 + k2] + d[26]*ddm6[234 + 3*k1 + k2] + da[k2]*dm3[k1] + da[k1]*dm3[k2] + db[k2]*dm4[k1] + db[k1]*dm4[k2] + db[3 + k2]*dm4[3 + k1] + db[3 + k1]*dm4[3 + k2] + db[6 + k2]*dm4[6 + k1] + db[6 + k1]*dm4[6 + k2] + dc[k2]*dm5[k1] + dc[k1]*dm5[k2] + 2*dc[3 + k2]*dm5[3 + k1] + 2*dc[3 + k1]*dm5[3 + k2] + 2*dc[6 + k2]*dm5[6 + k1] + 2*dc[6 + k1]*dm5[6 + k2] + dc[12 + k2]*dm5[12 + k1] + dc[12 + k1]*dm5[12 + k2] + 2*dc[15 + k2]*dm5[15 + k1] + 2*dc[15 + k1]*dm5[15 + k2] + dc[24 + k2]*dm5[24 + k1] + dc[24 + k1]*dm5[24 + k2] + dd[k2]*dm6[k1] + dd[k1]*dm6[k2] + 3*dd[3 + k2]*dm6[3 + k1] + 3*dd[3 + k1]*dm6[3 + k2] + 3*dd[6 + k2]*dm6[6 + k1] + 3*dd[6 + k1]*dm6[6 + k2] + 3*dd[12 + k2]*dm6[12 + k1] + 3*dd[12 + k1]*dm6[12 + k2] + 6*dd[15 + k2]*dm6[15 + k1] + 6*dd[15 + k1]*dm6[15 + k2] + 3*dd[24 + k2]*dm6[24 + k1] + 3*dd[24 + k1]*dm6[24 + k2] + dd[39 + k2]*dm6[39 + k1] + dd[39 + k1]*dm6[39 + k2] + 3*dd[42 + k2]*dm6[42 + k1] + 3*dd[42 + k1]*dm6[42 + k2] + 3*dd[51 + k2]*dm6[51 + k1] + 3*dd[51 + k1]*dm6[51 + k2] + dd[78 + k2]*dm6[78 + k1] + dd[78 + k1]*dm6[78 + k2]);
      rhs(11) = -(a*ddm3[9 + 3*k1 + k2] + b[0]*ddm4[9 + 3*k1 + k2] + b[1]*ddm4[36 + 3*k1 + k2] + b[2]*ddm4[45 + 3*k1 + k2] + c[0]*ddm5[9 + 3*k1 + k2] + 2*c[1]*ddm5[36 + 3*k1 + k2] + 2*c[2]*ddm5[45 + 3*k1 + k2] + c[4]*ddm5[117 + 3*k1 + k2] + 2*c[5]*ddm5[126 + 3*k1 + k2] + c[8]*ddm5[153 + 3*k1 + k2] + d[0]*ddm6[9 + 3*k1 + k2] + 3*d[1]*ddm6[36 + 3*k1 + k2] + 3*d[2]*ddm6[45 + 3*k1 + k2] + 3*d[4]*ddm6[117 + 3*k1 + k2] + 6*d[5]*ddm6[126 + 3*k1 + k2] + 3*d[8]*ddm6[153 + 3*k1 + k2] + d[13]*ddm6[360 + 3*k1 + k2] + 3*d[14]*ddm6[369 + 3*k1 + k2] + 3*d[17]*ddm6[396 + 3*k1 + k2] + d[26]*ddm6[477 + 3*k1 + k2] + da[k2]*dm3[3 + k1] + da[k1]*dm3[3 + k2] + db[k2]*dm4[3 + k1] + db[k1]*dm4[3 + k2] + db[3 + k2]*dm4[12 + k1] + db[3 + k1]*dm4[12 + k2] + db[6 + k2]*dm4[15 + k1] + db[6 + k1]*dm4[15 + k2] + dc[k2]*dm5[3 + k1] + dc[k1]*dm5[3 + k2] + 2*dc[3 + k2]*dm5[12 + k1] + 2*dc[3 + k1]*dm5[12 + k2] + 2*dc[6 + k2]*dm5[15 + k1] + 2*dc[6 + k1]*dm5[15 + k2] + dc[12 + k2]*dm5[39 + k1] + dc[12 + k1]*dm5[39 + k2] + 2*dc[15 + k2]*dm5[42 + k1] + 2*dc[15 + k1]*dm5[42 + k2] + dc[24 + k2]*dm5[51 + k1] + dc[24 + k1]*dm5[51 + k2] + dd[k2]*dm6[3 + k1] + dd[k1]*dm6[3 + k2] + 3*dd[3 + k2]*dm6[12 + k1] + 3*dd[3 + k1]*dm6[12 + k2] + 3*dd[6 + k2]*dm6[15 + k1] + 3*dd[6 + k1]*dm6[15 + k2] + 3*dd[12 + k2]*dm6[39 + k1] + 3*dd[12 + k1]*dm6[39 + k2] + 6*dd[15 + k2]*dm6[42 + k1] + 6*dd[15 + k1]*dm6[42 + k2] + 3*dd[24 + k2]*dm6[51 + k1] + 3*dd[24 + k1]*dm6[51 + k2] + dd[39 + k2]*dm6[120 + k1] + dd[39 + k1]*dm6[120 + k2] + 3*dd[42 + k2]*dm6[123 + k1] + 3*dd[42 + k1]*dm6[123 + k2] + 3*dd[51 + k2]*dm6[132 + k1] + 3*dd[51 + k1]*dm6[132 + k2] + dd[78 + k2]*dm6[159 + k1] + dd[78 + k1]*dm6[159 + k2]);
      rhs(12) = -(a*ddm3[18 + 3*k1 + k2] + b[0]*ddm4[18 + 3*k1 + k2] + b[1]*ddm4[45 + 3*k1 + k2] + b[2]*ddm4[72 + 3*k1 + k2] + c[0]*ddm5[18 + 3*k1 + k2] + 2*c[1]*ddm5[45 + 3*k1 + k2] + 2*c[2]*ddm5[72 + 3*k1 + k2] + c[4]*ddm5[126 + 3*k1 + k2] + 2*c[5]*ddm5[153 + 3*k1 + k2] + c[8]*ddm5[234 + 3*k1 + k2] + d[0]*ddm6[18 + 3*k1 + k2] + 3*d[1]*ddm6[45 + 3*k1 + k2] + 3*d[2]*ddm6[72 + 3*k1 + k2] + 3*d[4]*ddm6[126 + 3*k1 + k2] + 6*d[5]*ddm6[153 + 3*k1 + k2] + 3*d[8]*ddm6[234 + 3*k1 + k2] + d[13]*ddm6[369 + 3*k1 + k2] + 3*d[14]*ddm6[396 + 3*k1 + k2] + 3*d[17]*ddm6[477 + 3*k1 + k2] + d[26]*ddm6[720 + 3*k1 + k2] + da[k2]*dm3[6 + k1] + da[k1]*dm3[6 + k2] + db[k2]*dm4[6 + k1] + db[k1]*dm4[6 + k2] + db[3 + k2]*dm4[15 + k1] + db[3 + k1]*dm4[15 + k2] + db[6 + k2]*dm4[24 + k1] + db[6 + k1]*dm4[24 + k2] + dc[k2]*dm5[6 + k1] + dc[k1]*dm5[6 + k2] + 2*dc[3 + k2]*dm5[15 + k1] + 2*dc[3 + k1]*dm5[15 + k2] + 2*dc[6 + k2]*dm5[24 + k1] + 2*dc[6 + k1]*dm5[24 + k2] + dc[12 + k2]*dm5[42 + k1] + dc[12 + k1]*dm5[42 + k2] + 2*dc[15 + k2]*dm5[51 + k1] + 2*dc[15 + k1]*dm5[51 + k2] + dc[24 + k2]*dm5[78 + k1] + dc[24 + k1]*dm5[78 + k2] + dd[k2]*dm6[6 + k1] + dd[k1]*dm6[6 + k2] + 3*dd[3 + k2]*dm6[15 + k1] + 3*dd[3 + k1]*dm6[15 + k2] + 3*dd[6 + k2]*dm6[24 + k1] + 3*dd[6 + k1]*dm6[24 + k2] + 3*dd[12 + k2]*dm6[42 + k1] + 3*dd[12 + k1]*dm6[42 + k2] + 6*dd[15 + k2]*dm6[51 + k1] + 6*dd[15 + k1]*dm6[51 + k2] + 3*dd[24 + k2]*dm6[78 + k1] + 3*dd[24 + k1]*dm6[78 + k2] + dd[39 + k2]*dm6[123 + k1] + dd[39 + k1]*dm6[123 + k2] + 3*dd[42 + k2]*dm6[132 + k1] + 3*dd[42 + k1]*dm6[132 + k2] + 3*dd[51 + k2]*dm6[159 + k1] + 3*dd[51 + k1]*dm6[159 + k2] + dd[78 + k2]*dm6[240 + k1] + dd[78 + k1]*dm6[240 + k2]);
      rhs(13) = -(a*ddm3[36 + 3*k1 + k2] + b[0]*ddm4[36 + 3*k1 + k2] + b[1]*ddm4[117 + 3*k1 + k2] + b[2]*ddm4[126 + 3*k1 + k2] + c[0]*ddm5[36 + 3*k1 + k2] + 2*c[1]*ddm5[117 + 3*k1 + k2] + 2*c[2]*ddm5[126 + 3*k1 + k2] + c[4]*ddm5[360 + 3*k1 + k2] + 2*c[5]*ddm5[369 + 3*k1 + k2] + c[8]*ddm5[396 + 3*k1 + k2] + d[0]*ddm6[36 + 3*k1 + k2] + 3*d[1]*ddm6[117 + 3*k1 + k2] + 3*d[2]*ddm6[126 + 3*k1 + k2] + 3*d[4]*ddm6[360 + 3*k1 + k2] + 6*d[5]*ddm6[369 + 3*k1 + k2] + 3*d[8]*ddm6[396 + 3*k1 + k2] + d[13]*ddm6[1089 + 3*k1 + k2] + 3*d[14]*ddm6[1098 + 3*k1 + k2] + 3*d[17]*ddm6[1125 + 3*k1 + k2] + d[26]*ddm6[1206 + 3*k1 + k2] + da[k2]*dm3[12 + k1] + da[k1]*dm3[12 + k2] + db[k2]*dm4[12 + k1] + db[k1]*dm4[12 + k2] + db[3 + k2]*dm4[39 + k1] + db[3 + k1]*dm4[39 + k2] + db[6 + k2]*dm4[42 + k1] + db[6 + k1]*dm4[42 + k2] + dc[k2]*dm5[12 + k1] + dc[k1]*dm5[12 + k2] + 2*dc[3 + k2]*dm5[39 + k1] + 2*dc[3 + k1]*dm5[39 + k2] + 2*dc[6 + k2]*dm5[42 + k1] + 2*dc[6 + k1]*dm5[42 + k2] + dc[12 + k2]*dm5[120 + k1] + dc[12 + k1]*dm5[120 + k2] + 2*dc[15 + k2]*dm5[123 + k1] + 2*dc[15 + k1]*dm5[123 + k2] + dc[24 + k2]*dm5[132 + k1] + dc[24 + k1]*dm5[132 + k2] + dd[k2]*dm6[12 + k1] + dd[k1]*dm6[12 + k2] + 3*dd[3 + k2]*dm6[39 + k1] + 3*dd[3 + k1]*dm6[39 + k2] + 3*dd[6 + k2]*dm6[42 + k1] + 3*dd[6 + k1]*dm6[42 + k2] + 3*dd[12 + k2]*dm6[120 + k1] + 3*dd[12 + k1]*dm6[120 + k2] + 6*dd[15 + k2]*dm6[123 + k1] + 6*dd[15 + k1]*dm6[123 + k2] + 3*dd[24 + k2]*dm6[132 + k1] + 3*dd[24 + k1]*dm6[132 + k2] + dd[39 + k2]*dm6[363 + k1] + dd[39 + k1]*dm6[363 + k2] + 3*dd[42 + k2]*dm6[366 + k1] + 3*dd[42 + k1]*dm6[366 + k2] + 3*dd[51 + k2]*dm6[375 + k1] + 3*dd[51 + k1]*dm6[375 + k2] + dd[78 + k2]*dm6[402 + k1] + dd[78 + k1]*dm6[402 + k2]);
      rhs(14) = -(a*ddm3[45 + 3*k1 + k2] + b[0]*ddm4[45 + 3*k1 + k2] + b[1]*ddm4[126 + 3*k1 + k2] + b[2]*ddm4[153 + 3*k1 + k2] + c[0]*ddm5[45 + 3*k1 + k2] + 2*c[1]*ddm5[126 + 3*k1 + k2] + 2*c[2]*ddm5[153 + 3*k1 + k2] + c[4]*ddm5[369 + 3*k1 + k2] + 2*c[5]*ddm5[396 + 3*k1 + k2] + c[8]*ddm5[477 + 3*k1 + k2] + d[0]*ddm6[45 + 3*k1 + k2] + 3*d[1]*ddm6[126 + 3*k1 + k2] + 3*d[2]*ddm6[153 + 3*k1 + k2] + 3*d[4]*ddm6[369 + 3*k1 + k2] + 6*d[5]*ddm6[396 + 3*k1 + k2] + 3*d[8]*ddm6[477 + 3*k1 + k2] + d[13]*ddm6[1098 + 3*k1 + k2] + 3*d[14]*ddm6[1125 + 3*k1 + k2] + 3*d[17]*ddm6[1206 + 3*k1 + k2] + d[26]*ddm6[1449 + 3*k1 + k2] + da[k2]*dm3[15 + k1] + da[k1]*dm3[15 + k2] + db[k2]*dm4[15 + k1] + db[k1]*dm4[15 + k2] + db[3 + k2]*dm4[42 + k1] + db[3 + k1]*dm4[42 + k2] + db[6 + k2]*dm4[51 + k1] + db[6 + k1]*dm4[51 + k2] + dc[k2]*dm5[15 + k1] + dc[k1]*dm5[15 + k2] + 2*dc[3 + k2]*dm5[42 + k1] + 2*dc[3 + k1]*dm5[42 + k2] + 2*dc[6 + k2]*dm5[51 + k1] + 2*dc[6 + k1]*dm5[51 + k2] + dc[12 + k2]*dm5[123 + k1] + dc[12 + k1]*dm5[123 + k2] + 2*dc[15 + k2]*dm5[132 + k1] + 2*dc[15 + k1]*dm5[132 + k2] + dc[24 + k2]*dm5[159 + k1] + dc[24 + k1]*dm5[159 + k2] + dd[k2]*dm6[15 + k1] + dd[k1]*dm6[15 + k2] + 3*dd[3 + k2]*dm6[42 + k1] + 3*dd[3 + k1]*dm6[42 + k2] + 3*dd[6 + k2]*dm6[51 + k1] + 3*dd[6 + k1]*dm6[51 + k2] + 3*dd[12 + k2]*dm6[123 + k1] + 3*dd[12 + k1]*dm6[123 + k2] + 6*dd[15 + k2]*dm6[132 + k1] + 6*dd[15 + k1]*dm6[132 + k2] + 3*dd[24 + k2]*dm6[159 + k1] + 3*dd[24 + k1]*dm6[159 + k2] + dd[39 + k2]*dm6[366 + k1] + dd[39 + k1]*dm6[366 + k2] + 3*dd[42 + k2]*dm6[375 + k1] + 3*dd[42 + k1]*dm6[375 + k2] + 3*dd[51 + k2]*dm6[402 + k1] + 3*dd[51 + k1]*dm6[402 + k2] + dd[78 + k2]*dm6[483 + k1] + dd[78 + k1]*dm6[483 + k2]);
      rhs(15) = -(a*ddm3[72 + 3*k1 + k2] + b[0]*ddm4[72 + 3*k1 + k2] + b[1]*ddm4[153 + 3*k1 + k2] + b[2]*ddm4[234 + 3*k1 + k2] + c[0]*ddm5[72 + 3*k1 + k2] + 2*c[1]*ddm5[153 + 3*k1 + k2] + 2*c[2]*ddm5[234 + 3*k1 + k2] + c[4]*ddm5[396 + 3*k1 + k2] + 2*c[5]*ddm5[477 + 3*k1 + k2] + c[8]*ddm5[720 + 3*k1 + k2] + d[0]*ddm6[72 + 3*k1 + k2] + 3*d[1]*ddm6[153 + 3*k1 + k2] + 3*d[2]*ddm6[234 + 3*k1 + k2] + 3*d[4]*ddm6[396 + 3*k1 + k2] + 6*d[5]*ddm6[477 + 3*k1 + k2] + 3*d[8]*ddm6[720 + 3*k1 + k2] + d[13]*ddm6[1125 + 3*k1 + k2] + 3*d[14]*ddm6[1206 + 3*k1 + k2] + 3*d[17]*ddm6[1449 + 3*k1 + k2] + d[26]*ddm6[2178 + 3*k1 + k2] + da[k2]*dm3[24 + k1] + da[k1]*dm3[24 + k2] + db[k2]*dm4[24 + k1] + db[k1]*dm4[24 + k2] + db[3 + k2]*dm4[51 + k1] + db[3 + k1]*dm4[51 + k2] + db[6 + k2]*dm4[78 + k1] + db[6 + k1]*dm4[78 + k2] + dc[k2]*dm5[24 + k1] + dc[k1]*dm5[24 + k2] + 2*dc[3 + k2]*dm5[51 + k1] + 2*dc[3 + k1]*dm5[51 + k2] + 2*dc[6 + k2]*dm5[78 + k1] + 2*dc[6 + k1]*dm5[78 + k2] + dc[12 + k2]*dm5[132 + k1] + dc[12 + k1]*dm5[132 + k2] + 2*dc[15 + k2]*dm5[159 + k1] + 2*dc[15 + k1]*dm5[159 + k2] + dc[24 + k2]*dm5[240 + k1] + dc[24 + k1]*dm5[240 + k2] + dd[k2]*dm6[24 + k1] + dd[k1]*dm6[24 + k2] + 3*dd[3 + k2]*dm6[51 + k1] + 3*dd[3 + k1]*dm6[51 + k2] + 3*dd[6 + k2]*dm6[78 + k1] + 3*dd[6 + k1]*dm6[78 + k2] + 3*dd[12 + k2]*dm6[132 + k1] + 3*dd[12 + k1]*dm6[132 + k2] + 6*dd[15 + k2]*dm6[159 + k1] + 6*dd[15 + k1]*dm6[159 + k2] + 3*dd[24 + k2]*dm6[240 + k1] + 3*dd[24 + k1]*dm6[240 + k2] + dd[39 + k2]*dm6[375 + k1] + dd[39 + k1]*dm6[375 + k2] + 3*dd[42 + k2]*dm6[402 + k1] + 3*dd[42 + k1]*dm6[402 + k2] + 3*dd[51 + k2]*dm6[483 + k1] + 3*dd[51 + k1]*dm6[483 + k2] + dd[78 + k2]*dm6[726 + k1] + dd[78 + k1]*dm6[726 + k2]);
      rhs(16) = -(a*ddm3[117 + 3*k1 + k2] + b[0]*ddm4[117 + 3*k1 + k2] + b[1]*ddm4[360 + 3*k1 + k2] + b[2]*ddm4[369 + 3*k1 + k2] + c[0]*ddm5[117 + 3*k1 + k2] + 2*c[1]*ddm5[360 + 3*k1 + k2] + 2*c[2]*ddm5[369 + 3*k1 + k2] + c[4]*ddm5[1089 + 3*k1 + k2] + 2*c[5]*ddm5[1098 + 3*k1 + k2] + c[8]*ddm5[1125 + 3*k1 + k2] + d[0]*ddm6[117 + 3*k1 + k2] + 3*d[1]*ddm6[360 + 3*k1 + k2] + 3*d[2]*ddm6[369 + 3*k1 + k2] + 3*d[4]*ddm6[1089 + 3*k1 + k2] + 6*d[5]*ddm6[1098 + 3*k1 + k2] + 3*d[8]*ddm6[1125 + 3*k1 + k2] + d[13]*ddm6[3276 + 3*k1 + k2] + 3*d[14]*ddm6[3285 + 3*k1 + k2] + 3*d[17]*ddm6[3312 + 3*k1 + k2] + d[26]*ddm6[3393 + 3*k1 + k2] + da[k2]*dm3[39 + k1] + da[k1]*dm3[39 + k2] + db[k2]*dm4[39 + k1] + db[k1]*dm4[39 + k2] + db[3 + k2]*dm4[120 + k1] + db[3 + k1]*dm4[120 + k2] + db[6 + k2]*dm4[123 + k1] + db[6 + k1]*dm4[123 + k2] + dc[k2]*dm5[39 + k1] + dc[k1]*dm5[39 + k2] + 2*dc[3 + k2]*dm5[120 + k1] + 2*dc[3 + k1]*dm5[120 + k2] + 2*dc[6 + k2]*dm5[123 + k1] + 2*dc[6 + k1]*dm5[123 + k2] + dc[12 + k2]*dm5[363 + k1] + dc[12 + k1]*dm5[363 + k2] + 2*dc[15 + k2]*dm5[366 + k1] + 2*dc[15 + k1]*dm5[366 + k2] + dc[24 + k2]*dm5[375 + k1] + dc[24 + k1]*dm5[375 + k2] + dd[k2]*dm6[39 + k1] + dd[k1]*dm6[39 + k2] + 3*dd[3 + k2]*dm6[120 + k1] + 3*dd[3 + k1]*dm6[120 + k2] + 3*dd[6 + k2]*dm6[123 + k1] + 3*dd[6 + k1]*dm6[123 + k2] + 3*dd[12 + k2]*dm6[363 + k1] + 3*dd[12 + k1]*dm6[363 + k2] + 6*dd[15 + k2]*dm6[366 + k1] + 6*dd[15 + k1]*dm6[366 + k2] + 3*dd[24 + k2]*dm6[375 + k1] + 3*dd[24 + k1]*dm6[375 + k2] + dd[39 + k2]*dm6[1092 + k1] + dd[39 + k1]*dm6[1092 + k2] + 3*dd[42 + k2]*dm6[1095 + k1] + 3*dd[42 + k1]*dm6[1095 + k2] + 3*dd[51 + k2]*dm6[1104 + k1] + 3*dd[51 + k1]*dm6[1104 + k2] + dd[78 + k2]*dm6[1131 + k1] + dd[78 + k1]*dm6[1131 + k2]);
      rhs(17) = -(a*ddm3[126 + 3*k1 + k2] + b[0]*ddm4[126 + 3*k1 + k2] + b[1]*ddm4[369 + 3*k1 + k2] + b[2]*ddm4[396 + 3*k1 + k2] + c[0]*ddm5[126 + 3*k1 + k2] + 2*c[1]*ddm5[369 + 3*k1 + k2] + 2*c[2]*ddm5[396 + 3*k1 + k2] + c[4]*ddm5[1098 + 3*k1 + k2] + 2*c[5]*ddm5[1125 + 3*k1 + k2] + c[8]*ddm5[1206 + 3*k1 + k2] + d[0]*ddm6[126 + 3*k1 + k2] + 3*d[1]*ddm6[369 + 3*k1 + k2] + 3*d[2]*ddm6[396 + 3*k1 + k2] + 3*d[4]*ddm6[1098 + 3*k1 + k2] + 6*d[5]*ddm6[1125 + 3*k1 + k2] + 3*d[8]*ddm6[1206 + 3*k1 + k2] + d[13]*ddm6[3285 + 3*k1 + k2] + 3*d[14]*ddm6[3312 + 3*k1 + k2] + 3*d[17]*ddm6[3393 + 3*k1 + k2] + d[26]*ddm6[3636 + 3*k1 + k2] + da[k2]*dm3[42 + k1] + da[k1]*dm3[42 + k2] + db[k2]*dm4[42 + k1] + db[k1]*dm4[42 + k2] + db[3 + k2]*dm4[123 + k1] + db[3 + k1]*dm4[123 + k2] + db[6 + k2]*dm4[132 + k1] + db[6 + k1]*dm4[132 + k2] + dc[k2]*dm5[42 + k1] + dc[k1]*dm5[42 + k2] + 2*dc[3 + k2]*dm5[123 + k1] + 2*dc[3 + k1]*dm5[123 + k2] + 2*dc[6 + k2]*dm5[132 + k1] + 2*dc[6 + k1]*dm5[132 + k2] + dc[12 + k2]*dm5[366 + k1] + dc[12 + k1]*dm5[366 + k2] + 2*dc[15 + k2]*dm5[375 + k1] + 2*dc[15 + k1]*dm5[375 + k2] + dc[24 + k2]*dm5[402 + k1] + dc[24 + k1]*dm5[402 + k2] + dd[k2]*dm6[42 + k1] + dd[k1]*dm6[42 + k2] + 3*dd[3 + k2]*dm6[123 + k1] + 3*dd[3 + k1]*dm6[123 + k2] + 3*dd[6 + k2]*dm6[132 + k1] + 3*dd[6 + k1]*dm6[132 + k2] + 3*dd[12 + k2]*dm6[366 + k1] + 3*dd[12 + k1]*dm6[366 + k2] + 6*dd[15 + k2]*dm6[375 + k1] + 6*dd[15 + k1]*dm6[375 + k2] + 3*dd[24 + k2]*dm6[402 + k1] + 3*dd[24 + k1]*dm6[402 + k2] + dd[39 + k2]*dm6[1095 + k1] + dd[39 + k1]*dm6[1095 + k2] + 3*dd[42 + k2]*dm6[1104 + k1] + 3*dd[42 + k1]*dm6[1104 + k2] + 3*dd[51 + k2]*dm6[1131 + k1] + 3*dd[51 + k1]*dm6[1131 + k2] + dd[78 + k2]*dm6[1212 + k1] + dd[78 + k1]*dm6[1212 + k2]);
      rhs(18) = -(a*ddm3[153 + 3*k1 + k2] + b[0]*ddm4[153 + 3*k1 + k2] + b[1]*ddm4[396 + 3*k1 + k2] + b[2]*ddm4[477 + 3*k1 + k2] + c[0]*ddm5[153 + 3*k1 + k2] + 2*c[1]*ddm5[396 + 3*k1 + k2] + 2*c[2]*ddm5[477 + 3*k1 + k2] + c[4]*ddm5[1125 + 3*k1 + k2] + 2*c[5]*ddm5[1206 + 3*k1 + k2] + c[8]*ddm5[1449 + 3*k1 + k2] + d[0]*ddm6[153 + 3*k1 + k2] + 3*d[1]*ddm6[396 + 3*k1 + k2] + 3*d[2]*ddm6[477 + 3*k1 + k2] + 3*d[4]*ddm6[1125 + 3*k1 + k2] + 6*d[5]*ddm6[1206 + 3*k1 + k2] + 3*d[8]*ddm6[1449 + 3*k1 + k2] + d[13]*ddm6[3312 + 3*k1 + k2] + 3*d[14]*ddm6[3393 + 3*k1 + k2] + 3*d[17]*ddm6[3636 + 3*k1 + k2] + d[26]*ddm6[4365 + 3*k1 + k2] + da[k2]*dm3[51 + k1] + da[k1]*dm3[51 + k2] + db[k2]*dm4[51 + k1] + db[k1]*dm4[51 + k2] + db[3 + k2]*dm4[132 + k1] + db[3 + k1]*dm4[132 + k2] + db[6 + k2]*dm4[159 + k1] + db[6 + k1]*dm4[159 + k2] + dc[k2]*dm5[51 + k1] + dc[k1]*dm5[51 + k2] + 2*dc[3 + k2]*dm5[132 + k1] + 2*dc[3 + k1]*dm5[132 + k2] + 2*dc[6 + k2]*dm5[159 + k1] + 2*dc[6 + k1]*dm5[159 + k2] + dc[12 + k2]*dm5[375 + k1] + dc[12 + k1]*dm5[375 + k2] + 2*dc[15 + k2]*dm5[402 + k1] + 2*dc[15 + k1]*dm5[402 + k2] + dc[24 + k2]*dm5[483 + k1] + dc[24 + k1]*dm5[483 + k2] + dd[k2]*dm6[51 + k1] + dd[k1]*dm6[51 + k2] + 3*dd[3 + k2]*dm6[132 + k1] + 3*dd[3 + k1]*dm6[132 + k2] + 3*dd[6 + k2]*dm6[159 + k1] + 3*dd[6 + k1]*dm6[159 + k2] + 3*dd[12 + k2]*dm6[375 + k1] + 3*dd[12 + k1]*dm6[375 + k2] + 6*dd[15 + k2]*dm6[402 + k1] + 6*dd[15 + k1]*dm6[402 + k2] + 3*dd[24 + k2]*dm6[483 + k1] + 3*dd[24 + k1]*dm6[483 + k2] + dd[39 + k2]*dm6[1104 + k1] + dd[39 + k1]*dm6[1104 + k2] + 3*dd[42 + k2]*dm6[1131 + k1] + 3*dd[42 + k1]*dm6[1131 + k2] + 3*dd[51 + k2]*dm6[1212 + k1] + 3*dd[51 + k1]*dm6[1212 + k2] + dd[78 + k2]*dm6[1455 + k1] + dd[78 + k1]*dm6[1455 + k2]);
      rhs(19) = -(a*ddm3[234 + 3*k1 + k2] + b[0]*ddm4[234 + 3*k1 + k2] + b[1]*ddm4[477 + 3*k1 + k2] + b[2]*ddm4[720 + 3*k1 + k2] + c[0]*ddm5[234 + 3*k1 + k2] + 2*c[1]*ddm5[477 + 3*k1 + k2] + 2*c[2]*ddm5[720 + 3*k1 + k2] + c[4]*ddm5[1206 + 3*k1 + k2] + 2*c[5]*ddm5[1449 + 3*k1 + k2] + c[8]*ddm5[2178 + 3*k1 + k2] + d[0]*ddm6[234 + 3*k1 + k2] + 3*d[1]*ddm6[477 + 3*k1 + k2] + 3*d[2]*ddm6[720 + 3*k1 + k2] + 3*d[4]*ddm6[1206 + 3*k1 + k2] + 6*d[5]*ddm6[1449 + 3*k1 + k2] + 3*d[8]*ddm6[2178 + 3*k1 + k2] + d[13]*ddm6[3393 + 3*k1 + k2] + 3*d[14]*ddm6[3636 + 3*k1 + k2] + 3*d[17]*ddm6[4365 + 3*k1 + k2] + d[26]*ddm6[6552 + 3*k1 + k2] + da[k2]*dm3[78 + k1] + da[k1]*dm3[78 + k2] + db[k2]*dm4[78 + k1] + db[k1]*dm4[78 + k2] + db[3 + k2]*dm4[159 + k1] + db[3 + k1]*dm4[159 + k2] + db[6 + k2]*dm4[240 + k1] + db[6 + k1]*dm4[240 + k2] + dc[k2]*dm5[78 + k1] + dc[k1]*dm5[78 + k2] + 2*dc[3 + k2]*dm5[159 + k1] + 2*dc[3 + k1]*dm5[159 + k2] + 2*dc[6 + k2]*dm5[240 + k1] + 2*dc[6 + k1]*dm5[240 + k2] + dc[12 + k2]*dm5[402 + k1] + dc[12 + k1]*dm5[402 + k2] + 2*dc[15 + k2]*dm5[483 + k1] + 2*dc[15 + k1]*dm5[483 + k2] + dc[24 + k2]*dm5[726 + k1] + dc[24 + k1]*dm5[726 + k2] + dd[k2]*dm6[78 + k1] + dd[k1]*dm6[78 + k2] + 3*dd[3 + k2]*dm6[159 + k1] + 3*dd[3 + k1]*dm6[159 + k2] + 3*dd[6 + k2]*dm6[240 + k1] + 3*dd[6 + k1]*dm6[240 + k2] + 3*dd[12 + k2]*dm6[402 + k1] + 3*dd[12 + k1]*dm6[402 + k2] + 6*dd[15 + k2]*dm6[483 + k1] + 6*dd[15 + k1]*dm6[483 + k2] + 3*dd[24 + k2]*dm6[726 + k1] + 3*dd[24 + k1]*dm6[726 + k2] + dd[39 + k2]*dm6[1131 + k1] + dd[39 + k1]*dm6[1131 + k2] + 3*dd[42 + k2]*dm6[1212 + k1] + 3*dd[42 + k1]*dm6[1212 + k2] + 3*dd[51 + k2]*dm6[1455 + k1] + 3*dd[51 + k1]*dm6[1455 + k2] + dd[78 + k2]*dm6[2184 + k1] + dd[78 + k1]*dm6[2184 + k2]);
 
      lhs = solver.solve(rhs);

      dda[3*k1 + k2] = lhs(0);
      ddb[3*k1 + k2] = lhs(1);
      ddb[9 + 3*k1 + k2] = lhs(2);
      ddb[18 + 3*k1 + k2] = lhs(3);
      ddc[3*k1 + k2] = lhs(4);
      ddc[9 + 3*k1 + k2] = lhs(5);
      ddc[18 + 3*k1 + k2] = lhs(6);
      ddc[36 + 3*k1 + k2] = lhs(7);
      ddc[45 + 3*k1 + k2] = lhs(8);
      ddc[72 + 3*k1 + k2] = lhs(9);
      ddd[3*k1 + k2] = lhs(10);
      ddd[9 + 3*k1 + k2] = lhs(11);
      ddd[18 + 3*k1 + k2] = lhs(12);
      ddd[36 + 3*k1 + k2] = lhs(13);
      ddd[45 + 3*k1 + k2] = lhs(14);
      ddd[72 + 3*k1 + k2] = lhs(15);
      ddd[117 + 3*k1 + k2] = lhs(16);
      ddd[126 + 3*k1 + k2] = lhs(17);
      ddd[153 + 3*k1 + k2] = lhs(18);
      ddd[234 + 3*k1 + k2] = lhs(19);
    }
  }
  
  // Fill symmetries
  fillSymmetries<Dim<3>>(ddb);
  fillSymmetries<Dim<3>>(c, dc, ddc);
  fillSymmetries<Dim<3>>(d, dd, ddd);
}

} // end namespace anonymous

  //------------------------------------------------------------------------------
  // Templated version for faster summation of moments
  //------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
void
computeRKCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                     const TableKernel<Dimension>& W,
                     const FieldList<Dimension, typename Dimension::Scalar>& volume,
                     const FieldList<Dimension, typename Dimension::Vector>& position,
                     const FieldList<Dimension, typename Dimension::SymTensor>& H,
                     FieldList<Dimension, typename Dimension::Scalar>& A,
                     FieldList<Dimension, typename Dimension::Vector>& B,
                     FieldList<Dimension, typename Dimension::Tensor>& C,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& D,
                     FieldList<Dimension, typename Dimension::Vector>& gradA,
                     FieldList<Dimension, typename Dimension::Tensor>& gradB,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& gradD,
                     FieldList<Dimension, typename Dimension::Tensor>& hessA,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& hessB,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& hessC,
                     FieldList<Dimension, typename Dimension::FifthRankTensor>& hessD) {
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;
  
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

  // Get an identity tensor
  Tensor identityTensor;
  identityTensor.Identity();

  // Add stuff to the moments for a single point combination (so self contribution is not duplicated)
  const bool etaInterp = false;
  auto addPointToMoments = [&](const int nodeListi, const int nodei,
                               const int nodeListj, const int nodej,
                               RKMomentValues<Dimension>& mom) {
    // Get H and eta values
    const auto xi = position(nodeListi, nodei);
    const auto xj = position(nodeListj, nodej);
    const auto xij = xi - xj;
    const auto Hij = H(nodeListj, nodej);
    const auto etaij = Hij * xij;
    const auto etaMagInvij = safeInv(etaij.magnitude());
    const auto H2ij = Hij.square();
    const auto Hetaij = Hij * etaij * etaMagInvij;
    const auto Heta2ij = Hetaij.selfdyad();
          
    // Get kernel values
    const auto Kij = W(etaij, Hij);
    const auto gradKij = W.grad(etaij, Hij);
    const auto hessKij = W.grad2(etaij, Hij);
    const auto Wij = Kij;
    const auto gradWij = Hij * etaij * etaMagInvij * gradKij;
    const auto hessWij = Tensor((H2ij - Heta2ij) * etaMagInvij * gradKij + Heta2ij * hessKij);
    const auto vj = volume(nodeListj, nodej);

    // Get function to reproduce and its derivative
    const auto g = etaInterp ? etaij : xij;
    const auto dg = etaInterp ? Tensor(Hij) : identityTensor;

    // Add the values to the moments
    addToMoments<Dimension, correctionOrder>(g, dg, Wij, gradWij, hessWij, vj, mom);
          
    return;
  };
  
  // Get the moments
  RKMomentValues<Dimension> rkMoments;

  // Compute things point by point
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = connectivityMap.numNodes(nodeListi);
    
    for (auto nodei = 0; nodei < numNodes; ++nodei) {
      
      // Zero out the moments
      rkMoments.template zero<correctionOrder>();
      
      // Calculate M values
      addPointToMoments(nodeListi, nodei, nodeListi, nodei, rkMoments); // self contribution
      const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, nodei);
      for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
        for (auto nodej : connectivity[nodeListj]) {
          addPointToMoments(nodeListi, nodei, nodeListj, nodej, rkMoments); // other points
        }
      }

      // Compute corrections
      switch (correctionOrder) {
      case CRKOrder::ZerothOrder:
        computeCorrections(rkMoments,
                           A(nodeListi, nodei), gradA(nodeListi, nodei), hessA(nodeListi, nodei));
        break;
      case CRKOrder::LinearOrder:
        computeCorrections(rkMoments,
                           A(nodeListi, nodei), B(nodeListi, nodei),
                           gradA(nodeListi, nodei), gradB(nodeListi, nodei),
                           hessA(nodeListi, nodei), hessB(nodeListi, nodei));
        break;
      case CRKOrder::QuadraticOrder:
        computeCorrections(rkMoments,
                           A(nodeListi, nodei), B(nodeListi, nodei), C(nodeListi, nodei),
                           gradA(nodeListi, nodei), gradB(nodeListi, nodei), gradC(nodeListi, nodei),
                           hessA(nodeListi, nodei), hessB(nodeListi, nodei), hessC(nodeListi, nodei));
        break;
      case CRKOrder::CubicOrder:
        computeCorrections(rkMoments,
                           A(nodeListi, nodei), B(nodeListi, nodei), C(nodeListi, nodei), D(nodeListi, nodei),
                           gradA(nodeListi, nodei), gradB(nodeListi, nodei), gradC(nodeListi, nodei), gradD(nodeListi, nodei),
                           hessA(nodeListi, nodei), hessB(nodeListi, nodei), hessC(nodeListi, nodei), hessD(nodeListi, nodei));
        break;
      }
    } // nodei
  } // nodeListi
} // computeRKCorrections

  //------------------------------------------------------------------------------
  // Choose correct templated version
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
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& D,
                     FieldList<Dimension, typename Dimension::Vector>& gradA,
                     FieldList<Dimension, typename Dimension::Tensor>& gradB,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& gradC,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& gradD,
                     FieldList<Dimension, typename Dimension::Tensor>& hessA,
                     FieldList<Dimension, typename Dimension::ThirdRankTensor>& hessB,
                     FieldList<Dimension, typename Dimension::FourthRankTensor>& hessC,
                     FieldList<Dimension, typename Dimension::FifthRankTensor>& hessD) {
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    return computeRKCorrections<Dimension, CRKOrder::ZerothOrder>(connectivityMap,
                                                                  W, volume, position, H,
                                                                  A, B, C, D,
                                                                  gradA, gradB, gradC, gradD,
                                                                  hessA, hessB, hessC, hessD);
  case CRKOrder::LinearOrder:
    return computeRKCorrections<Dimension, CRKOrder::LinearOrder>(connectivityMap,
                                                                  W, volume, position, H,
                                                                  A, B, C, D,
                                                                  gradA, gradB, gradC, gradD,
                                                                  hessA, hessB, hessC, hessD);
  case CRKOrder::QuadraticOrder:
    return computeRKCorrections<Dimension, CRKOrder::QuadraticOrder>(connectivityMap,
                                                                     W, volume, position, H,
                                                                     A, B, C, D,
                                                                     gradA, gradB, gradC, gradD,
                                                                     hessA, hessB, hessC, hessD);
  case CRKOrder::CubicOrder:
    return computeRKCorrections<Dimension, CRKOrder::CubicOrder>(connectivityMap,
                                                                 W, volume, position, H,
                                                                 A, B, C, D,
                                                                 gradA, gradB, gradC, gradD,
                                                                 hessA, hessB, hessC, hessD);
  }
}

// //------------------------------------------------------------------------------
// // Do instantiations explicitly for now
// //------------------------------------------------------------------------------
// template<>
// void
// computeRKCorrections(const ConnectivityMap<Dim<1>>& connectivityMap,
//                      const TableKernel<Dim<1>>& W,
//                      const FieldList<Dim<1>, Dim<1>::Scalar>& volume,
//                      const FieldList<Dim<1>, Dim<1>::Vector>& position,
//                      const FieldList<Dim<1>, Dim<1>::SymTensor>& H,
//                      const CRKOrder correctionOrder,
//                      FieldList<Dim<1>, Dim<1>::Scalar>& A,
//                      FieldList<Dim<1>, Dim<1>::Vector>& B,
//                      FieldList<Dim<1>, Dim<1>::Tensor>& C,
//                      FieldList<Dim<1>, Dim<1>::ThirdRankTensor>& D,
//                      FieldList<Dim<1>, Dim<1>::Vector>& gradA,
//                      FieldList<Dim<1>, Dim<1>::Tensor>& gradB,
//                      FieldList<Dim<1>, Dim<1>::ThirdRankTensor>& gradC,
//                      FieldList<Dim<1>, Dim<1>::FourthRankTensor>& gradD,
//                      FieldList<Dim<1>, Dim<1>::Tensor>& hessA,
//                      FieldList<Dim<1>, Dim<1>::ThirdRankTensor>& hessB,
//                      FieldList<Dim<1>, Dim<1>::FourthRankTensor>& hessC,
//                      FieldList<Dim<1>, Dim<1>::FifthRankTensor>& hessD);

// template<>
// void
// computeRKCorrections(const ConnectivityMap<Dim<2>>& connectivityMap,
//                      const TableKernel<Dim<2>>& W,
//                      const FieldList<Dim<2>, Dim<2>::Scalar>& volume,
//                      const FieldList<Dim<2>, Dim<2>::Vector>& position,
//                      const FieldList<Dim<2>, Dim<2>::SymTensor>& H,
//                      const CRKOrder correctionOrder,
//                      FieldList<Dim<2>, Dim<2>::Scalar>& A,
//                      FieldList<Dim<2>, Dim<2>::Vector>& B,
//                      FieldList<Dim<2>, Dim<2>::Tensor>& C,
//                      FieldList<Dim<2>, Dim<2>::ThirdRankTensor>& D,
//                      FieldList<Dim<2>, Dim<2>::Vector>& gradA,
//                      FieldList<Dim<2>, Dim<2>::Tensor>& gradB,
//                      FieldList<Dim<2>, Dim<2>::ThirdRankTensor>& gradC,
//                      FieldList<Dim<2>, Dim<2>::FourthRankTensor>& gradD,
//                      FieldList<Dim<2>, Dim<2>::Tensor>& hessA,
//                      FieldList<Dim<2>, Dim<2>::ThirdRankTensor>& hessB,
//                      FieldList<Dim<2>, Dim<2>::FourthRankTensor>& hessC,
//                      FieldList<Dim<2>, Dim<2>::FifthRankTensor>& hessD);

// template<>
// void
// computeRKCorrections(const ConnectivityMap<Dim<3>>& connectivityMap,
//                      const TableKernel<Dim<3>>& W,
//                      const FieldList<Dim<3>, Dim<3>::Scalar>& volume,
//                      const FieldList<Dim<3>, Dim<3>::Vector>& position,
//                      const FieldList<Dim<3>, Dim<3>::SymTensor>& H,
//                      const CRKOrder correctionOrder,
//                      FieldList<Dim<3>, Dim<3>::Scalar>& A,
//                      FieldList<Dim<3>, Dim<3>::Vector>& B,
//                      FieldList<Dim<3>, Dim<3>::Tensor>& C,
//                      FieldList<Dim<3>, Dim<3>::ThirdRankTensor>& D,
//                      FieldList<Dim<3>, Dim<3>::Vector>& gradA,
//                      FieldList<Dim<3>, Dim<3>::Tensor>& gradB,
//                      FieldList<Dim<3>, Dim<3>::ThirdRankTensor>& gradC,
//                      FieldList<Dim<3>, Dim<3>::FourthRankTensor>& gradD,
//                      FieldList<Dim<3>, Dim<3>::Tensor>& hessA,
//                      FieldList<Dim<3>, Dim<3>::ThirdRankTensor>& hessB,
//                      FieldList<Dim<3>, Dim<3>::FourthRankTensor>& hessC,
//                      FieldList<Dim<3>, Dim<3>::FifthRankTensor>& hessD);

} // end namespace Spheral
