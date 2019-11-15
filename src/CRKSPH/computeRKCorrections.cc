//------------------------------------------------------------------------------
// Compute the RK correction terms
//------------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include "computeRKCorrections.hh"
#include "Eigen/Dense"
#include "Field/Field.hh"
#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "NodeList/NodeList.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Helper class for holding the moments
//------------------------------------------------------------------------------
template<typename Dimension>
struct RKMomentValues {
  RKMomentValues()
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
  void zero() {
    std::fill(std::begin(m0), std::end(m0), 0);
    std::fill(std::begin(m1), std::end(m1), 0);
    std::fill(std::begin(m2), std::end(m2), 0);
    std::fill(std::begin(m3), std::end(m3), 0);
    std::fill(std::begin(m4), std::end(m4), 0);
    std::fill(std::begin(m5), std::end(m5), 0);
    std::fill(std::begin(m6), std::end(m6), 0);
    std::fill(std::begin(dm0), std::end(dm0), 0);
    std::fill(std::begin(dm1), std::end(dm1), 0);
    std::fill(std::begin(dm2), std::end(dm2), 0);
    std::fill(std::begin(dm3), std::end(dm3), 0);
    std::fill(std::begin(dm4), std::end(dm4), 0);
    std::fill(std::begin(dm5), std::end(dm5), 0);
    std::fill(std::begin(dm6), std::end(dm6), 0);
    std::fill(std::begin(ddm0), std::end(ddm0), 0);
    std::fill(std::begin(ddm1), std::end(ddm1), 0);
    std::fill(std::begin(ddm2), std::end(ddm2), 0);
    std::fill(std::begin(ddm3), std::end(ddm3), 0);
    std::fill(std::begin(ddm4), std::end(ddm4), 0);
    std::fill(std::begin(ddm5), std::end(ddm5), 0);
    std::fill(std::begin(ddm6), std::end(ddm6), 0);
  }
};

//------------------------------------------------------------------------------
// Add value to moments
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
void
addToMoments(const typename Dimension::Vector& eta,
             const typename Dimension::Tensor& deta,
             const typename Dimension::ThirdRankTensor& ddeta,
             const typename Dimension::Scalar& w,
             const typename Dimension::Vector& dw,
             const typename Dimension::Tensor& ddw,
             const typename Dimension::Scalar& v,
             RKMomentValues<Dimension>& moments);

// d=1, o=0
template<>
void
addToMoments<Dim<1>, CRKOrder::ZerothOrder>(const Dim<1>::Vector& eta,
                                            const Dim<1>::Tensor& deta,
                                            const Dim<1>::ThirdRankTensor& ddeta,
                                            const Dim<1>::Scalar& w,
                                            const Dim<1>::Vector& dw,
                                            const Dim<1>::Tensor& ddw,
                                            const Dim<1>::Scalar& v,
                                            RKMomentValues<Dim<1>>& moments) {
  const auto dim = Dim<1>::nDim;
  auto& m0 = moments.m0;
  auto& dm0 = moments.dm0;
  auto& ddm0 = moments.ddm0;
  
  // Moments
  m0 += v*w;
  // Gradients
  dm0[0] += v*dw[0];
  // Hessians
  ddm0[0] += v*ddw[0];
}
// d=2, o=0
template<>
void
addToMoments<Dim<2>, CRKOrder::ZerothOrder>(const Dim<2>::Vector& eta,
                                            const Dim<2>::Tensor& deta,
                                            const Dim<2>::ThirdRankTensor& ddeta,
                                            const Dim<2>::Scalar& w,
                                            const Dim<2>::Vector& dw,
                                            const Dim<2>::Tensor& ddw,
                                            const Dim<2>::Scalar& v,
                                            RKMomentValues<Dim<2>>& moments) {
  const auto dim = Dim<2>::nDim;
  auto& m0 = moments.m0;
  auto& dm0 = moments.dm0;
  auto& ddm0 = moments.ddm0;
  
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
void
addToMoments<Dim<3>, CRKOrder::ZerothOrder>(const Dim<3>::Vector& eta,
                                            const Dim<3>::Tensor& deta,
                                            const Dim<3>::ThirdRankTensor& ddeta,
                                            const Dim<3>::Scalar& w,
                                            const Dim<3>::Vector& dw,
                                            const Dim<3>::Tensor& ddw,
                                            const Dim<3>::Scalar& v,
                                            RKMomentValues<Dim<3>>& moments) {
  const auto dim = Dim<3>::nDim;
  auto& m0 = moments.m0;
  auto& dm0 = moments.dm0;
  auto& ddm0 = moments.ddm0;
  
  // Moments
  m0 += v*w;
  // Gradients
  for (auto k1 = 0; k1 < dim; ++k1) {
    dm0[k1] += v*dw[k1];
  }
  // Hessians
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      ddm0[3*k1 + k2] += v*ddw[3*k1 + k2]
    }
  }
}

//------------------------------------------------------------------------------
// Compute the corrections for a single point
//------------------------------------------------------------------------------
template<typename Dimension>
void
computeCorrections(const RKMomentValues<Dimension>& moments,
                   typename Dimension::Scalar& a,
                   typename Dimension::Vector& da,
                   typename Dimension::Tensor& dda);
template<typename Dimension>
void
computeCorrections(const RKMomentValues<Dimension>& moments,
                   typename Dimension::Scalar& a,
                   typename Dimension::Vector& b,
                   typename Dimension::Vector& da,
                   typename Dimension::Tensor& db,
                   typename Dimension::Tensor& dda,
                   typename Dimension::ThirdRankTensor& ddb);
template<typename Dimension>
void
computeCorrections(const RKMomentValues<Dimension>& moments,
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
void
computeCorrections(const RKMomentValues<Dimension>& moments,
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
                   typename Dimension::FifthRankTensor& ddd);
// d=1, o=0
template<>
void
computeCorrections(const RKMomentValues<<Dim<1>>& moments,
                   Dim<1>::Scalar& a,
                   Dim<1>::Vector& da,
                   Dim<1>::Tensor& dda) {
  const auto dim = Dim<1>::nDim;
  const auto& m0 = moments.m0;
  const auto& dm0 = moments.dm0;
  const auto& ddm0 = moments.ddm0;
  
  // Value
  a[0] = 1./m0;
  // Gradients
  da[0] = -a*dm0[k1]/m0;
  // Hessians
  dda[0] = -(a*ddm0[k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2])/m0;
}

// d=2, o=0
template<>
void
computeCorrections(const RKMomentValues<<Dim<1>>& moments,
                   Dim<2>::Scalar& a,
                   Dim<2>::Vector& da,
                   Dim<2>::Tensor& dda) {
  
  const auto dim = Dim<2>::nDim;
  const auto& m0 = moments.m0;
  const auto& dm0 = moments.dm0;
  const auto& ddm0 = moments.ddm0;
  
  // Value
  a[0] = 1./m0;
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
void
computeCorrections(const RKMomentValues<<Dim<1>>& moments,
                   Dim<2>::Scalar& a,
                   Dim<2>::Vector& da,
                   Dim<2>::Tensor& dda) {
  
  const auto dim = Dim<3>::nDim;
  const auto& m0 = moments.m0;
  const auto& dm0 = moments.dm0;
  const auto& ddm0 = moments.ddm0;
  
  // Value
  a[0] = 1./m0;
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
void
computeCorrections(const RKMomentValues<Dimension>& moments,
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
  const auto& m0 = moments.m0;
  const auto& m1 = moments.m1;
  const auto& m2 = moments.m2;
  const auto& dm0 = moments.dm0;
  const auto& dm1 = moments.dm1;
  const auto& dm2 = moments.dm2;
  const auto& ddm0 = moments.ddm0;
  const auto& ddm1 = moments.ddm1;
  const auto& ddm2 = moments.ddm2;
  const auto k1 = 0;
  const auto k2 = 0;
  const auto p1 = 0;

  // Get matrix to invert
  VectorType rhs;
  MatrixType matrix;
  VectorType lhs;
  matrix(0, 0) = m0;
  matrix(1, 0) = m1[0];
  matrix(0, 1) = m1[0];
  matrix(1, 1) = m2[0];
  auto solver = matrix.colPivHouseholderQr();

  // Solve for values
  rhs(0) = 1.;
  rhs(1) = 0.;
  lhs = solver.solve(rhs);
  a[0] = lhs(0);
  b[0] = lhs(1);
  // Solve for derivatives
  rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1]);
  rhs(1) = -(a*dm1[k1 + p1] + b[0]*dm2[k1 + p1]);
  lhs = solver.solve(rhs);
  da[0] = lhs(0);
  db[0] = lhs(1);
  // Solve for second derivatives
  rhs(0) = -(a*ddm0[k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2]);
  rhs(1) = -(a*ddm1[k1 + k2 + p1] + b[0]*ddm2[k1 + k2 + p1] + da[k2]*dm1[k1 + p1] + da[k1]*dm1[k2 + p1] + db[k2]*dm2[k1 + p1] + db[k1]*dm2[k2 + p1]);
  lhs = solver.solve(rhs);
  da[0] = lhs(0);
  db[0] = lhs(1);
}
// d=2, o=1
template<>
void
computeCorrections(const RKMomentValues<Dimension>& moments,
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
  const auto& m0 = moments.m0;
  const auto& m1 = moments.m1;
  const auto& m2 = moments.m2;
  const auto& dm0 = moments.dm0;
  const auto& dm1 = moments.dm1;
  const auto& dm2 = moments.dm2;
  const auto& ddm0 = moments.ddm0;
  const auto& ddm1 = moments.ddm1;
  const auto& ddm2 = moments.ddm2;

  // Get matrix to invert
  VectorType rhs;
  MatrixType matrix;
  VectorType lhs;
  matrix(0, 0) = m0;
  for (auto q1 = 0; q1 < dim; ++q1) {
    matrix(0, q1+1) = m1[q1];
    matrix(q1+1, 0) = m1[q1];
    for (auto q2 = 0; q2 < dim; ++q2) {
      matrix(q1+1,q2+1) = m2[q2*2 + q1];
    }
  }
  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = 1.;
  for (auto p1 = 0; p1 < dim; ++p1) {
    const auto row = p1 + 1;
    rhs(row) = 0.;
  }
  lhs = solver.solve(rhs);
  a[0] = lhs(0);
  for (auto p1 = 0; p1 < dim; ++p1) {
    const auto row = p1 + 1;
    b[p1] = lhs(row);
  }
  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[2 + k1]);
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      rhs(row) = -(a*dm1[k1 + 2*p1] + b[0]*dm2[k1 + 2*p1] + b[1]*dm2[4 + k1 + 2*p1]); 
    }
    lhs = solver.solve(rhs);
    da[k1] = rhs(0);
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      db[2*p1 + k1] = rhs(row);
    }
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[2*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[2*k1 + k2] + b[1]*ddm1[4 + 2*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[2 + k2]*dm1[2 + k1] + db[2 + k1]*dm1[2 + k2]);
      for (auto p1 = 0; p1 < dim; ++p1) {
        const auto row = p1 + 1;
        rhs(row) = -(a*ddm1[2*k1 + k2 + 4*p1] + b[0]*ddm2[2*k1 + k2 + 4*p1] + b[1]*ddm2[8 + 2*k1 + k2 + 4*p1] + da[k2]*dm1[k1 + 2*p1] + da[k1]*dm1[k2 + 2*p1] + db[k2]*dm2[k1 + 2*p1] + db[k1]*dm2[k2 + 2*p1] + db[2 + k2]*dm2[4 + k1 + 2*p1] + db[2 + k1]*dm2[4 + k2 + 2*p1]);
      }
      lhs = solver.solve(rhs);
      dda[2*k1 + k2] = solution(0);
      for (auto p1 = 0; p1 < dim; ++p1) {
        const auto row = p1 + 1;
        ddb[4*p1 + 2*k1 + k2] = solution(p1);
      }
    }
  }
}
// d=3, o=1
template<>
void
computeCorrections(const RKMomentValues<Dimension>& moments,
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
  const auto& m0 = moments.m0;
  const auto& m1 = moments.m1;
  const auto& m2 = moments.m2;
  const auto& dm0 = moments.dm0;
  const auto& dm1 = moments.dm1;
  const auto& dm2 = moments.dm2;
  const auto& ddm0 = moments.ddm0;
  const auto& ddm1 = moments.ddm1;
  const auto& ddm2 = moments.ddm2;

  // Get matrix to invert
  VectorType rhs;
  MatrixType matrix;
  VectorType lhs;
  matrix(0, 0) = m0;
  for (auto q1 = 0; q1 < dim; ++q1) {
    matrix(0, q1+1) = m1[q1];
    matrix(q1+1, 0) = m1[q1];
    for (auto q2 = 0; q2 < dim; ++q2) {
      matrix(q1+1,q2+1) = m2[q2*3 + q1];
    }
  }
  auto solver = matrix.colPivHouseholderQr();
  
  // Solve for values
  rhs(0) = 1.;
  for (auto p1 = 0; p1 < dim; ++p1) {
    const auto row = p1 + 1;
    rhs(row) = 0.;
  }
  lhs = solver.solve(rhs);
  a[0] = lhs(0);
  for (auto p1 = 0; p1 < dim; ++p1) {
    const auto row = p1 + 1;
    b[p1] = lhs(row);
  }
  // Solve for derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    rhs(0) = -(a*dm0[k1] + b[0]*dm1[k1] + b[1]*dm1[3 + k1] + b[2]*dm1[6 + k1]);
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      rhs(row) = -(a*dm1[k1 + 3*p1] + b[0]*dm2[k1 + 3*p1] + b[1]*dm2[9 + k1 + 3*p1] + b[2]*dm2[18 + k1 + 3*p1]);
    }
    lhs = solver.solve(rhs);
    da[k1] = rhs(0);
    for (auto p1 = 0; p1 < dim; ++p1) {
      const auto row = p1 + 1;
      db[3*p1 + k1] = rhs(row);
    }
  }
  // Solve for second derivatives
  for (auto k1 = 0; k1 < dim; ++k1) {
    for (auto k2 = k1; k2 < dim; ++k2) {
      rhs(0) = -(a*ddm0[3*k1 + k2] + da[k2]*dm0[k1] + da[k1]*dm0[k2] + b[0]*ddm1[3*k1 + k2] + b[1]*ddm1[9 + 3*k1 + k2] + b[2]*ddm1[18 + 3*k1 + k2] + db[k2]*dm1[k1] + db[k1]*dm1[k2] + db[3 + k2]*dm1[3 + k1] + db[3 + k1]*dm1[3 + k2] + db[6 + k2]*dm1[6 + k1] + db[6 + k1]*dm1[6 + k2]);
      for (auto p1 = 0; p1 < dim; ++p1) {
        const auto row = p1 + 1;
        rhs(row) = -(a*ddm1[3*k1 + k2 + 9*p1] + b[0]*ddm2[3*k1 + k2 + 9*p1] + b[1]*ddm2[27 + 3*k1 + k2 + 9*p1] + b[2]*ddm2[54 + 3*k1 + k2 + 9*p1] + da[k2]*dm1[k1 + 3*p1] + da[k1]*dm1[k2 + 3*p1] + db[k2]*dm2[k1 + 3*p1] + db[k1]*dm2[k2 + 3*p1] + db[3 + k2]*dm2[9 + k1 + 3*p1] + db[3 + k1]*dm2[9 + k2 + 3*p1] + db[6 + k2]*dm2[18 + k1 + 3*p1] + db[6 + k1]*dm2[18 + k2 + 3*p1]);
      }
      lhs = solver.solve(rhs);
      dda[3*k1 + k2] = solution(0);
      for (auto p1 = 0; p1 < dim; ++p1) {
        const auto row = p1 + 1;
        ddb[6*p1 + 3*k1 + k2] = solution(p1);
      }
    }
  }
}

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
      moments.zero();
      
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

          addToMoments<Dimension, correctionOrder>(eta, gradEta, hessEta, Wij, gradWij, hessWij, vj, moments);
        }
      }
      
      // Compute corrections
      switch (correctionOrder) {
      case CRKOrder::ZerothOrder:
        computeCorrections(moments,
                           A(nodeListi, nodei), gradA(nodeListi, nodei), hessA(nodeListi, nodei));
        break;
      case CRKOrder::LinearOrder:
        computeCorrections(moments,
                           A(nodeListi, nodei), B(nodeListi, nodei),
                           gradA(nodeListi, nodei), gradB(nodeListi, nodei),
                           hessA(nodeListi, nodei), hessB(nodeListi, nodei));
          break;
      case CRKOrder::QuadraticOrder:
        computeCorrections(moments,
                           A(nodeListi, nodei), B(nodeListi, nodei), C(nodeListi, nodei),
                           gradA(nodeListi, nodei), gradB(nodeListi, nodei), gradC(nodeListi, nodei),
                           hessA(nodeListi, nodei), hessB(nodeListi, nodei), hessC(nodeListi, nodei));
        break;
      case CRKOrder::CubicOrder:
        computeCorrections(moments,
                           A(nodeListi, nodei), B(nodeListi, nodei), C(nodeListi, nodei), D(nodeListi, nodei),
                           gradA(nodeListi, nodei), gradB(nodeListi, nodei), gradC(nodeListi, nodei), gradD(nodeListi, nodei),
                           hessA(nodeListi, nodei), hessB(nodeListi, nodei), hessC(nodeListi, nodei), hessD(nodeListi, nodei));
        break;
      }
    } // nodei
  } // nodeListi
} // computeRKCorrections

//------------------------------------------------------------------------------
// Call the correct templated version
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
  // case CRKOrder::QuadraticOrder:
  //   return computeRKCorrections<Dimension, CRKOrder::QuadraticOrder>(connectivityMap,
  //                                                                 W, volume, position, H,
  //                                                                 A, B, C, D,
  //                                                                 gradA, gradB, gradC, gradD,
  //                                                                 hessA, hessB, hessC, hessD);
  // case CRKOrder::CubicOrder:
  //   return computeRKCorrections<Dimension, CRKOrder::CubicOrder>(connectivityMap,
  //                                                                 W, volume, position, H,
  //                                                                 A, B, C, D,
  //                                                                 gradA, gradB, gradC, gradD,
  //                                                                 hessA, hessB, hessC, hessD);
  default:
    ASSERT2(false, "order not implemented");
  }
}

} // end namespace spheral
