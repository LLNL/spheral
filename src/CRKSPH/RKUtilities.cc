//---------------------------------Spheral++----------------------------------//
// RKUtilities
//
// Evaluate the RK kernel
//----------------------------------------------------------------------------//
#include "RKUtilities.hh"

#include "Kernel/TableKernel.hh"
#include "CRKSPHCorrectionParams.hh"
#include "Geometry/Dimension.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the corrected kernel value
//------------------------------------------------------------------------------
template<>
typename Dim<1>::Scalar
evaluateRKKernel<Dim<1>>(const TableKernel<Dim<1>>& kernel,
                         const CRKOrder correctionOrder,
                         const Dim<1>::Vector& eta,
                         const Dim<1>::SymTensor& H,
                         const Dim<1>::Vector& g,
                         const Dim<1>::Scalar& a,
                         const Dim<1>::Vector& b,
                         const Dim<1>::Tensor& c,
                         const Dim<1>::ThirdRankTensor& d) {
  // Evaluate the kernel
  const auto w = kernel(eta, H);
  
  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    return a*w;
  case CRKOrder::LinearOrder:
    return w*(a + b[0]*g[0]);
  case CRKOrder::QuadraticOrder:
    return w*(a + g[0]*(b[0] + c[0]*g[0]));
  case CRKOrder::CubicOrder:
    return w*(a + g[0]*(b[0] + g[0]*(c[0] + d[0]*g[0])));
  }
}

template<>
typename Dim<2>::Scalar
evaluateRKKernel<Dim<2>>(const TableKernel<Dim<2>>& kernel,
                         const CRKOrder correctionOrder,
                         const Dim<2>::Vector& eta,
                         const Dim<2>::SymTensor& H,
                         const Dim<2>::Vector& g,
                         const Dim<2>::Scalar& a,
                         const Dim<2>::Vector& b,
                         const Dim<2>::Tensor& c,
                         const Dim<2>::ThirdRankTensor& d) {
  // Evaluate the kernel
  const auto w = kernel(eta, H);
  
  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    return a*w;
  case CRKOrder::LinearOrder:
    return w*(a + b[0]*g[0] + b[1]*g[1]);
  case CRKOrder::QuadraticOrder:
    return w*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + g[1]*(b[1] + (c[1] + c[2])*g[0] + c[3]*g[1]));
  case CRKOrder::CubicOrder:
    return w*(a + g[0]*(b[0] + g[0]*(c[0] + d[0]*g[0])) + (b[1] + g[0]*(c[1] + c[2] + (d[1] + d[2] + d[4])*g[0]))*g[1] + (c[3] + (d[3] + d[5] + d[6])*g[0])*std::pow(g[1],2) + d[7]*std::pow(g[1],3));
  }
}

template<>
typename Dim<3>::Scalar
evaluateRKKernel<Dim<3>>(const TableKernel<Dim<3>>& kernel,
                         const CRKOrder correctionOrder,
                         const Dim<3>::Vector& eta,
                         const Dim<3>::SymTensor& H,
                         const Dim<3>::Vector& g,
                         const Dim<3>::Scalar& a,
                         const Dim<3>::Vector& b,
                         const Dim<3>::Tensor& c,
                         const Dim<3>::ThirdRankTensor& d) {
  // Evaluate the kernel
  const auto w = kernel(eta, H);
  
  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    return a*w;
  case CRKOrder::LinearOrder:
    return w*(a + b[0]*g[0] + b[1]*g[1] + b[2]*g[2]);
  case CRKOrder::QuadraticOrder:
    return w*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + g[1]*(b[1] + (c[1] + c[3])*g[0] + c[4]*g[1]) + (b[2] + (c[2] + c[6])*g[0] + (c[5] + c[7])*g[1])*g[2] + c[8]*std::pow(g[2],2));
  case CRKOrder::CubicOrder:
    return w*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + d[0]*std::pow(g[0],3) + b[1]*g[1] + c[1]*g[0]*g[1] + c[3]*g[0]*g[1] + d[1]*std::pow(g[0],2)*g[1] + d[3]*std::pow(g[0],2)*g[1] + d[9]*std::pow(g[0],2)*g[1] + c[4]*std::pow(g[1],2) + d[4]*g[0]*std::pow(g[1],2) + d[10]*g[0]*std::pow(g[1],2) + d[12]*g[0]*std::pow(g[1],2) + d[13]*std::pow(g[1],3) + (b[2] + g[0]*(c[2] + c[6] + (d[2] + d[6] + d[18])*g[0]) + (c[5] + c[7] + (d[5] + d[7] + d[11] + d[15] + d[19] + d[21])*g[0])*g[1] + (d[14] + d[16] + d[22])*std::pow(g[1],2))*g[2] + (c[8] + (d[8] + d[20] + d[24])*g[0] + (d[17] + d[23] + d[25])*g[1])*std::pow(g[2],2) + d[26]*std::pow(g[2],3));
  }
}

//------------------------------------------------------------------------------
// Compute the corrected kernel gradient
//------------------------------------------------------------------------------
template<>
Dim<1>::Vector
evaluateRKGradient<Dim<1>>(const TableKernel<Dim<1>>& kernel,
                           const CRKOrder correctionOrder,
                           const Dim<1>::Vector& eta,
                           const Dim<1>::SymTensor& H,
                           const Dim<1>::Vector& g,
                           const Dim<1>::Tensor& dg,
                           const Dim<1>::Scalar& a,
                           const Dim<1>::Vector& b,
                           const Dim<1>::Tensor& c,
                           const Dim<1>::ThirdRankTensor& d,
                           const Dim<1>::Vector& da,
                           const Dim<1>::Tensor& db,
                           const Dim<1>::ThirdRankTensor& dc,
                           const Dim<1>::FourthRankTensor& dd) {
  typedef Dim<1>::Vector Vector;
  typedef Dim<1>::Tensor Tensor;
  
  // Evaluate the kernel
  const auto deta = Tensor(H);
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto k = kernel(eta, H);
  const auto dk = kernel.grad(eta, H);
  const auto w = k;
  const auto dw = Heta * dk;
  
  // Initialize values
  Vector grad = Vector::zero;
  const auto dim = Dim<1>::nDim;
  const auto k1 = 0;
  
  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    grad[k1] = w*da[k1] + a*dw[k1];
    break;
  case CRKOrder::LinearOrder:
    grad[k1] = dw[k1]*(a + b[0]*g[0]) + w*(da[k1] + g[0]*db[k1] + b[0]*dg[k1]);
    break;
  case CRKOrder::QuadraticOrder:
    grad[k1] = dw[k1]*(a + g[0]*(b[0] + c[0]*g[0])) + w*(da[k1] + g[0]*(db[k1] + g[0]*dc[k1]) + (b[0] + 2*c[0]*g[0])*dg[k1]);
    break;
  case CRKOrder::CubicOrder:
    grad[k1] = dw[k1]*(a + g[0]*(b[0] + g[0]*(c[0] + d[0]*g[0]))) + w*(da[k1] + b[0]*dg[k1] + g[0]*(db[k1] + 2*c[0]*dg[k1] + g[0]*(dc[k1] + g[0]*dd[k1] + 3*d[0]*dg[k1])));
    break;
  }
  return grad;
}

template<>
Dim<2>::Vector
evaluateRKGradient<Dim<2>>(const TableKernel<Dim<2>>& kernel,
                           const CRKOrder correctionOrder,
                           const Dim<2>::Vector& eta,
                           const Dim<2>::SymTensor& H,
                           const Dim<2>::Vector& g,
                           const Dim<2>::Tensor& dg,
                           const Dim<2>::Scalar& a,
                           const Dim<2>::Vector& b,
                           const Dim<2>::Tensor& c,
                           const Dim<2>::ThirdRankTensor& d,
                           const Dim<2>::Vector& da,
                           const Dim<2>::Tensor& db,
                           const Dim<2>::ThirdRankTensor& dc,
                           const Dim<2>::FourthRankTensor& dd) {
  typedef Dim<2>::Vector Vector;
  typedef Dim<2>::Tensor Tensor;
  
  // Evaluate the kernel
  const auto deta = Tensor(H);
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto k = kernel(eta, H);
  const auto dk = kernel.grad(eta, H);
  const auto w = k;
  const auto dw = Heta * dk;
  
  // Initialize values
  Vector grad = Vector::zero;
  const auto dim = Dim<2>::nDim;
  
  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = w*da[k1] + a*dw[k1];
    }
    break;
  case CRKOrder::LinearOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = dw[k1]*(a + b[0]*g[0] + b[1]*g[1]) + w*(da[k1] + g[0]*db[k1] + g[1]*db[2 + k1] + b[0]*dg[k1] + b[1]*dg[2 + k1]);
    }
    break;
  case CRKOrder::QuadraticOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = dw[k1]*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + g[1]*(b[1] + (c[1] + c[2])*g[0] + c[3]*g[1])) + w*(da[k1] + g[0]*db[k1] + g[1]*db[2 + k1] + std::pow(g[0],2)*dc[k1] + g[0]*g[1]*dc[2 + k1] + g[0]*g[1]*dc[4 + k1] + std::pow(g[1],2)*dc[6 + k1] + b[0]*dg[k1] + 2*c[0]*g[0]*dg[k1] + c[1]*g[1]*dg[k1] + c[2]*g[1]*dg[k1] + (b[1] + (c[1] + c[2])*g[0] + 2*c[3]*g[1])*dg[2 + k1]);
    }
    break;
  case CRKOrder::CubicOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = dw[k1]*(a + g[0]*(b[0] + g[0]*(c[0] + d[0]*g[0])) + (b[1] + g[0]*(c[1] + c[2] + (d[1] + d[2] + d[4])*g[0]))*g[1] + (c[3] + (d[3] + d[5] + d[6])*g[0])*std::pow(g[1],2) + d[7]*std::pow(g[1],3)) + w*(da[k1] + g[0]*db[k1] + g[1]*db[2 + k1] + std::pow(g[0],2)*dc[k1] + g[0]*g[1]*dc[2 + k1] + g[0]*g[1]*dc[4 + k1] + std::pow(g[1],2)*dc[6 + k1] + std::pow(g[0],3)*dd[k1] + std::pow(g[0],2)*g[1]*dd[2 + k1] + std::pow(g[0],2)*g[1]*dd[4 + k1] + g[0]*std::pow(g[1],2)*dd[6 + k1] + std::pow(g[0],2)*g[1]*dd[8 + k1] + g[0]*std::pow(g[1],2)*dd[10 + k1] + g[0]*std::pow(g[1],2)*dd[12 + k1] + std::pow(g[1],3)*dd[14 + k1] + b[0]*dg[k1] + 2*c[0]*g[0]*dg[k1] + 3*d[0]*std::pow(g[0],2)*dg[k1] + c[1]*g[1]*dg[k1] + c[2]*g[1]*dg[k1] + 2*d[1]*g[0]*g[1]*dg[k1] + 2*d[2]*g[0]*g[1]*dg[k1] + 2*d[4]*g[0]*g[1]*dg[k1] + d[3]*std::pow(g[1],2)*dg[k1] + d[5]*std::pow(g[1],2)*dg[k1] + d[6]*std::pow(g[1],2)*dg[k1] + (b[1] + g[0]*(c[1] + c[2] + (d[1] + d[2] + d[4])*g[0]) + 2*(c[3] + (d[3] + d[5] + d[6])*g[0])*g[1] + 3*d[7]*std::pow(g[1],2))*dg[2 + k1]);
    }
    break;
  }
  return grad;
}

template<>
Dim<3>::Vector
evaluateRKGradient<Dim<3>>(const TableKernel<Dim<3>>& kernel,
                           const CRKOrder correctionOrder,
                           const Dim<3>::Vector& eta,
                           const Dim<3>::SymTensor& H,
                           const Dim<3>::Vector& g,
                           const Dim<3>::Tensor& dg,
                           const Dim<3>::Scalar& a,
                           const Dim<3>::Vector& b,
                           const Dim<3>::Tensor& c,
                           const Dim<3>::ThirdRankTensor& d,
                           const Dim<3>::Vector& da,
                           const Dim<3>::Tensor& db,
                           const Dim<3>::ThirdRankTensor& dc,
                           const Dim<3>::FourthRankTensor& dd) {
  typedef Dim<3>::Vector Vector;
  typedef Dim<3>::Tensor Tensor;
  
  // Evaluate the kernel
  const auto deta = Tensor(H);
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto k = kernel(eta, H);
  const auto dk = kernel.grad(eta, H);
  const auto w = k;
  const auto dw = Heta * dk;
  
  // Initialize values
  Vector grad = Vector::zero;
  const auto dim = Dim<3>::nDim;
  
  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = w*da[k1] + a*dw[k1];
    }
    break;
  case CRKOrder::LinearOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = w*(a + b[0]*g[0] + b[1]*g[1] + b[2]*g[2]);
    }
    break;
  case CRKOrder::QuadraticOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = dw[k1]*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + g[1]*(b[1] + (c[1] + c[3])*g[0] + c[4]*g[1]) + (b[2] + (c[2] + c[6])*g[0] + (c[5] + c[7])*g[1])*g[2] + c[8]*std::pow(g[2],2)) + w*(da[k1] + g[0]*db[k1] + g[1]*db[3 + k1] + g[2]*db[6 + k1] + std::pow(g[0],2)*dc[k1] + g[0]*g[1]*dc[3 + k1] + g[0]*g[2]*dc[6 + k1] + g[0]*g[1]*dc[9 + k1] + std::pow(g[1],2)*dc[12 + k1] + g[1]*g[2]*dc[15 + k1] + g[0]*g[2]*dc[18 + k1] + g[1]*g[2]*dc[21 + k1] + std::pow(g[2],2)*dc[24 + k1] + b[0]*dg[k1] + 2*c[0]*g[0]*dg[k1] + c[1]*g[1]*dg[k1] + c[3]*g[1]*dg[k1] + c[2]*g[2]*dg[k1] + c[6]*g[2]*dg[k1] + b[1]*dg[3 + k1] + c[1]*g[0]*dg[3 + k1] + c[3]*g[0]*dg[3 + k1] + 2*c[4]*g[1]*dg[3 + k1] + c[5]*g[2]*dg[3 + k1] + c[7]*g[2]*dg[3 + k1] + (b[2] + (c[2] + c[6])*g[0] + (c[5] + c[7])*g[1] + 2*c[8]*g[2])*dg[6 + k1]);
    }
    break;
  case CRKOrder::CubicOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] = dw[k1]*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + d[0]*std::pow(g[0],3) + b[1]*g[1] + c[1]*g[0]*g[1] + c[3]*g[0]*g[1] + d[1]*std::pow(g[0],2)*g[1] + d[3]*std::pow(g[0],2)*g[1] + d[9]*std::pow(g[0],2)*g[1] + c[4]*std::pow(g[1],2) + d[4]*g[0]*std::pow(g[1],2) + d[10]*g[0]*std::pow(g[1],2) + d[12]*g[0]*std::pow(g[1],2) + d[13]*std::pow(g[1],3) + (b[2] + g[0]*(c[2] + c[6] + (d[2] + d[6] + d[18])*g[0]) + (c[5] + c[7] + (d[5] + d[7] + d[11] + d[15] + d[19] + d[21])*g[0])*g[1] + (d[14] + d[16] + d[22])*std::pow(g[1],2))*g[2] + (c[8] + (d[8] + d[20] + d[24])*g[0] + (d[17] + d[23] + d[25])*g[1])*std::pow(g[2],2) + d[26]*std::pow(g[2],3)) + w*(da[k1] + g[0]*db[k1] + g[1]*db[3 + k1] + g[2]*db[6 + k1] + std::pow(g[0],2)*dc[k1] + g[0]*g[1]*dc[3 + k1] + g[0]*g[2]*dc[6 + k1] + g[0]*g[1]*dc[9 + k1] + std::pow(g[1],2)*dc[12 + k1] + g[1]*g[2]*dc[15 + k1] + g[0]*g[2]*dc[18 + k1] + g[1]*g[2]*dc[21 + k1] + std::pow(g[2],2)*dc[24 + k1] + std::pow(g[0],3)*dd[k1] + std::pow(g[0],2)*g[1]*dd[3 + k1] + std::pow(g[0],2)*g[2]*dd[6 + k1] + std::pow(g[0],2)*g[1]*dd[9 + k1] + g[0]*std::pow(g[1],2)*dd[12 + k1] + g[0]*g[1]*g[2]*dd[15 + k1] + std::pow(g[0],2)*g[2]*dd[18 + k1] + g[0]*g[1]*g[2]*dd[21 + k1] + g[0]*std::pow(g[2],2)*dd[24 + k1] + std::pow(g[0],2)*g[1]*dd[27 + k1] + g[0]*std::pow(g[1],2)*dd[30 + k1] + g[0]*g[1]*g[2]*dd[33 + k1] + g[0]*std::pow(g[1],2)*dd[36 + k1] + std::pow(g[1],3)*dd[39 + k1] + std::pow(g[1],2)*g[2]*dd[42 + k1] + g[0]*g[1]*g[2]*dd[45 + k1] + std::pow(g[1],2)*g[2]*dd[48 + k1] + g[1]*std::pow(g[2],2)*dd[51 + k1] + std::pow(g[0],2)*g[2]*dd[54 + k1] + g[0]*g[1]*g[2]*dd[57 + k1] + g[0]*std::pow(g[2],2)*dd[60 + k1] + g[0]*g[1]*g[2]*dd[63 + k1] + std::pow(g[1],2)*g[2]*dd[66 + k1] + g[1]*std::pow(g[2],2)*dd[69 + k1] + g[0]*std::pow(g[2],2)*dd[72 + k1] + g[1]*std::pow(g[2],2)*dd[75 + k1] + std::pow(g[2],3)*dd[78 + k1] + b[0]*dg[k1] + 2*c[0]*g[0]*dg[k1] + 3*d[0]*std::pow(g[0],2)*dg[k1] + c[1]*g[1]*dg[k1] + c[3]*g[1]*dg[k1] + 2*d[1]*g[0]*g[1]*dg[k1] + 2*d[3]*g[0]*g[1]*dg[k1] + 2*d[9]*g[0]*g[1]*dg[k1] + d[4]*std::pow(g[1],2)*dg[k1] + d[10]*std::pow(g[1],2)*dg[k1] + d[12]*std::pow(g[1],2)*dg[k1] + c[2]*g[2]*dg[k1] + c[6]*g[2]*dg[k1] + 2*d[2]*g[0]*g[2]*dg[k1] + 2*d[6]*g[0]*g[2]*dg[k1] + 2*d[18]*g[0]*g[2]*dg[k1] + d[5]*g[1]*g[2]*dg[k1] + d[7]*g[1]*g[2]*dg[k1] + d[11]*g[1]*g[2]*dg[k1] + d[15]*g[1]*g[2]*dg[k1] + d[19]*g[1]*g[2]*dg[k1] + d[21]*g[1]*g[2]*dg[k1] + d[8]*std::pow(g[2],2)*dg[k1] + d[20]*std::pow(g[2],2)*dg[k1] + d[24]*std::pow(g[2],2)*dg[k1] + b[1]*dg[3 + k1] + c[1]*g[0]*dg[3 + k1] + c[3]*g[0]*dg[3 + k1] + d[1]*std::pow(g[0],2)*dg[3 + k1] + d[3]*std::pow(g[0],2)*dg[3 + k1] + d[9]*std::pow(g[0],2)*dg[3 + k1] + 2*c[4]*g[1]*dg[3 + k1] + 2*d[4]*g[0]*g[1]*dg[3 + k1] + 2*d[10]*g[0]*g[1]*dg[3 + k1] + 2*d[12]*g[0]*g[1]*dg[3 + k1] + 3*d[13]*std::pow(g[1],2)*dg[3 + k1] + c[5]*g[2]*dg[3 + k1] + c[7]*g[2]*dg[3 + k1] + d[5]*g[0]*g[2]*dg[3 + k1] + d[7]*g[0]*g[2]*dg[3 + k1] + d[11]*g[0]*g[2]*dg[3 + k1] + d[15]*g[0]*g[2]*dg[3 + k1] + d[19]*g[0]*g[2]*dg[3 + k1] + d[21]*g[0]*g[2]*dg[3 + k1] + 2*d[14]*g[1]*g[2]*dg[3 + k1] + 2*d[16]*g[1]*g[2]*dg[3 + k1] + 2*d[22]*g[1]*g[2]*dg[3 + k1] + d[17]*std::pow(g[2],2)*dg[3 + k1] + d[23]*std::pow(g[2],2)*dg[3 + k1] + d[25]*std::pow(g[2],2)*dg[3 + k1] + (b[2] + c[2]*g[0] + c[6]*g[0] + d[2]*std::pow(g[0],2) + d[6]*std::pow(g[0],2) + d[18]*std::pow(g[0],2) + c[5]*g[1] + c[7]*g[1] + d[5]*g[0]*g[1] + d[7]*g[0]*g[1] + d[11]*g[0]*g[1] + d[15]*g[0]*g[1] + d[19]*g[0]*g[1] + d[21]*g[0]*g[1] + d[14]*std::pow(g[1],2) + d[16]*std::pow(g[1],2) + d[22]*std::pow(g[1],2) + 2*(c[8] + (d[8] + d[20] + d[24])*g[0] + (d[17] + d[23] + d[25])*g[1])*g[2] + 3*d[26]*std::pow(g[2],2))*dg[6 + k1]);
    }
    break;
  }
  return grad;
}

//------------------------------------------------------------------------------
// Compute the corrected kernel gradient
//------------------------------------------------------------------------------
template<>
Dim<1>::SymTensor
evaluateRKHessian<Dim<1>>(const TableKernel<Dim<1>>& kernel,
                          const CRKOrder correctionOrder,
                          const Dim<1>::Vector& eta,
                          const Dim<1>::SymTensor& H,
                          const Dim<1>::Vector& g,
                          const Dim<1>::Tensor& dg,
                          const Dim<1>::Scalar& a,
                          const Dim<1>::Vector& b,
                          const Dim<1>::Tensor& c,
                          const Dim<1>::ThirdRankTensor& d,
                          const Dim<1>::Vector& da,
                          const Dim<1>::Tensor& db,
                          const Dim<1>::ThirdRankTensor& dc,
                          const Dim<1>::FourthRankTensor& dd,
                          const Dim<1>::Tensor& dda,
                          const Dim<1>::ThirdRankTensor& ddb,
                          const Dim<1>::FourthRankTensor& ddc,
                          const Dim<1>::FifthRankTensor& ddd) {
  typedef Dim<1>::Tensor Tensor;
  typedef Dim<1>::SymTensor SymTensor;
  typedef Dim<1>::ThirdRankTensor ThirdRankTensor;

  // Evaluate the kernel
  const auto deta = Tensor(H);
  const auto ddeta = ThirdRankTensor::zero;
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto k = kernel(eta, H);
  const auto dk = kernel.grad(eta, H);
  const auto ddk = kernel.grad2(eta, H);
  const auto w = k;
  const auto dw = Heta * dk;
  const auto ddw = Tensor((H2 - Heta2) * etaMagInv * dk + Heta2 * ddk);
  
  // Initialize values
  Tensor hess = Tensor::zero;
  const auto dim = Dim<1>::nDim;
  const auto k1 = 0;
  const auto k2 = 0;

  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    hess[k1 + k2] = w*dda[k1 + k2] + a*ddw[k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2];
    break;
  case CRKOrder::LinearOrder:
    hess[k1 + k2] = ddw[k1 + k2]*(a + b[0]*g[0]) + dw[k2]*(da[k1] + g[0]*db[k1] + b[0]*dg[k1]) + dw[k1]*(da[k2] + g[0]*db[k2] + b[0]*dg[k2]) + w*(dda[k1 + k2] + g[0]*ddb[k1 + k2] + db[k2]*dg[k1] + db[k1]*dg[k2]);
    break;
  case CRKOrder::QuadraticOrder:
    hess[k1 + k2]  = w*dda[k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2] + ddw[k1 + k2]*(a + g[0]*(b[0] + c[0]*g[0])) + w*g[0]*ddb[k1 + k2] + w*std::pow(g[0],2)*ddc[k1 + k2] + 2*c[0]*dw[k2]*g[0]*dg[k1] + w*db[k2]*dg[k1] + 2*w*g[0]*dc[k2]*dg[k1] + dw[k2]*(g[0]*(db[k1] + g[0]*dc[k1]) + b[0]*dg[k1]) + (2*c[0]*dw[k1]*g[0] + w*(db[k1] + 2*g[0]*dc[k1] + 2*c[0]*dg[k1]))*dg[k2] + dw[k1]*(g[0]*(db[k2] + g[0]*dc[k2]) + b[0]*dg[k2]);
    break;
  case CRKOrder::CubicOrder:
    hess[k1 + k2] = w*dda[k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2] + ddw[k1 + k2]*(a + g[0]*(b[0] + g[0]*(c[0] + d[0]*g[0]))) + w*g[0]*ddb[k1 + k2] + w*std::pow(g[0],2)*ddc[k1 + k2] + w*std::pow(g[0],3)*ddd[k1 + k2] + 2*c[0]*dw[k2]*g[0]*dg[k1] + 3*d[0]*dw[k2]*std::pow(g[0],2)*dg[k1] + w*db[k2]*dg[k1] + 2*w*g[0]*dc[k2]*dg[k1] + 3*w*std::pow(g[0],2)*dd[k2]*dg[k1] + dw[k2]*(g[0]*(db[k1] + g[0]*(dc[k1] + g[0]*dd[k1])) + b[0]*dg[k1]) + (dw[k1]*g[0]*(2*c[0] + 3*d[0]*g[0]) + w*(db[k1] + g[0]*(2*dc[k1] + 3*g[0]*dd[k1]) + 2*(c[0] + 3*d[0]*g[0])*dg[k1]))*dg[k2] + dw[k1]*(g[0]*(db[k2] + g[0]*(dc[k2] + g[0]*dd[k2])) + b[0]*dg[k2]);
    break;
  }
  return SymTensor(hess);
}

template<>
Dim<2>::SymTensor
evaluateRKHessian<Dim<2>>(const TableKernel<Dim<2>>& kernel,
                          const CRKOrder correctionOrder,
                          const Dim<2>::Vector& eta,
                          const Dim<2>::SymTensor& H,
                          const Dim<2>::Vector& g,
                          const Dim<2>::Tensor& dg,
                          const Dim<2>::Scalar& a,
                          const Dim<2>::Vector& b,
                          const Dim<2>::Tensor& c,
                          const Dim<2>::ThirdRankTensor& d,
                          const Dim<2>::Vector& da,
                          const Dim<2>::Tensor& db,
                          const Dim<2>::ThirdRankTensor& dc,
                          const Dim<2>::FourthRankTensor& dd,
                          const Dim<2>::Tensor& dda,
                          const Dim<2>::ThirdRankTensor& ddb,
                          const Dim<2>::FourthRankTensor& ddc,
                          const Dim<2>::FifthRankTensor& ddd) {
  typedef Dim<2>::Tensor Tensor;
  typedef Dim<2>::SymTensor SymTensor;
  typedef Dim<2>::ThirdRankTensor ThirdRankTensor;

  // Evaluate the kernel
  const auto deta = Tensor(H);
  const auto ddeta = ThirdRankTensor::zero;
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto k = kernel(eta, H);
  const auto dk = kernel.grad(eta, H);
  const auto ddk = kernel.grad2(eta, H);
  const auto w = k;
  const auto dw = Heta * dk;
  const auto ddw = Tensor((H2 - Heta2) * etaMagInv * dk + Heta2 * ddk);
  
  // Initialize values
  Tensor hess = Tensor::zero;
  const auto dim = Dim<2>::nDim;

  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[2*k1 + k2] = w*dda[2*k1 + k2] + a*ddw[2*k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2];
      }
    }
    break;
  case CRKOrder::LinearOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[2*k1 + k2] = ddw[2*k1 + k2]*(a + b[0]*g[0] + b[1]*g[1]) + dw[k2]*(da[k1] + g[0]*db[k1] + g[1]*db[2 + k1] + b[0]*dg[k1] + b[1]*dg[2 + k1]) + dw[k1]*(da[k2] + g[0]*db[k2] + g[1]*db[2 + k2] + b[0]*dg[k2] + b[1]*dg[2 + k2]) + w*(dda[2*k1 + k2] + g[0]*ddb[2*k1 + k2] + g[1]*ddb[4 + 2*k1 + k2] + db[k2]*dg[k1] + db[k1]*dg[k2] + db[2 + k2]*dg[2 + k1] + db[2 + k1]*dg[2 + k2]);
      }
    }
    break;
  case CRKOrder::QuadraticOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[2*k1 + k2] = w*dda[2*k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2] + ddw[2*k1 + k2]*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + g[1]*(b[1] + (c[1] + c[2])*g[0] + c[3]*g[1])) + w*g[0]*ddb[2*k1 + k2] + w*g[1]*ddb[4 + 2*k1 + k2] + w*std::pow(g[0],2)*ddc[2*k1 + k2] + w*g[0]*g[1]*ddc[4 + 2*k1 + k2] + w*g[0]*g[1]*ddc[8 + 2*k1 + k2] + w*std::pow(g[1],2)*ddc[12 + 2*k1 + k2] + 2*c[0]*dw[k2]*g[0]*dg[k1] + c[1]*dw[k2]*g[1]*dg[k1] + c[2]*dw[k2]*g[1]*dg[k1] + w*db[k2]*dg[k1] + 2*w*g[0]*dc[k2]*dg[k1] + w*g[1]*dc[2 + k2]*dg[k1] + w*g[1]*dc[4 + k2]*dg[k1] + dw[k2]*(g[0]*db[k1] + std::pow(g[0],2)*dc[k1] + g[1]*(db[2 + k1] + g[0]*(dc[2 + k1] + dc[4 + k1]) + g[1]*dc[6 + k1]) + b[0]*dg[k1]) + 2*c[0]*dw[k1]*g[0]*dg[k2] + c[1]*dw[k1]*g[1]*dg[k2] + c[2]*dw[k1]*g[1]*dg[k2] + w*db[k1]*dg[k2] + 2*w*g[0]*dc[k1]*dg[k2] + w*g[1]*dc[2 + k1]*dg[k2] + w*g[1]*dc[4 + k1]*dg[k2] + 2*w*c[0]*dg[k1]*dg[k2] + dw[k1]*(g[0]*db[k2] + std::pow(g[0],2)*dc[k2] + g[1]*(db[2 + k2] + g[0]*(dc[2 + k2] + dc[4 + k2]) + g[1]*dc[6 + k2]) + b[0]*dg[k2]) + b[1]*dw[k2]*dg[2 + k1] + c[1]*dw[k2]*g[0]*dg[2 + k1] + c[2]*dw[k2]*g[0]*dg[2 + k1] + 2*c[3]*dw[k2]*g[1]*dg[2 + k1] + w*db[2 + k2]*dg[2 + k1] + w*g[0]*dc[2 + k2]*dg[2 + k1] + w*g[0]*dc[4 + k2]*dg[2 + k1] + 2*w*g[1]*dc[6 + k2]*dg[2 + k1] + w*c[1]*dg[k2]*dg[2 + k1] + w*c[2]*dg[k2]*dg[2 + k1] + (dw[k1]*(b[1] + (c[1] + c[2])*g[0] + 2*c[3]*g[1]) + w*(db[2 + k1] + g[0]*(dc[2 + k1] + dc[4 + k1]) + 2*g[1]*dc[6 + k1] + (c[1] + c[2])*dg[k1] + 2*c[3]*dg[2 + k1]))*dg[2 + k2];
      }
    }
    break;
  case CRKOrder::CubicOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[2*k1 + k2] = w*dda[2*k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2] + ddw[2*k1 + k2]*(a + g[0]*(b[0] + g[0]*(c[0] + d[0]*g[0])) + (b[1] + g[0]*(c[1] + c[2] + (d[1] + d[2] + d[4])*g[0]))*g[1] + (c[3] + (d[3] + d[5] + d[6])*g[0])*std::pow(g[1],2) + d[7]*std::pow(g[1],3)) + w*g[0]*ddb[2*k1 + k2] + w*g[1]*ddb[4 + 2*k1 + k2] + w*std::pow(g[0],2)*ddc[2*k1 + k2] + w*g[0]*g[1]*ddc[4 + 2*k1 + k2] + w*g[0]*g[1]*ddc[8 + 2*k1 + k2] + w*std::pow(g[1],2)*ddc[12 + 2*k1 + k2] + w*std::pow(g[0],3)*ddd[2*k1 + k2] + w*std::pow(g[0],2)*g[1]*ddd[4 + 2*k1 + k2] + w*std::pow(g[0],2)*g[1]*ddd[8 + 2*k1 + k2] + w*g[0]*std::pow(g[1],2)*ddd[12 + 2*k1 + k2] + w*std::pow(g[0],2)*g[1]*ddd[16 + 2*k1 + k2] + w*g[0]*std::pow(g[1],2)*ddd[20 + 2*k1 + k2] + w*g[0]*std::pow(g[1],2)*ddd[24 + 2*k1 + k2] + w*std::pow(g[1],3)*ddd[28 + 2*k1 + k2] + 2*c[0]*dw[k2]*g[0]*dg[k1] + 3*d[0]*dw[k2]*std::pow(g[0],2)*dg[k1] + c[1]*dw[k2]*g[1]*dg[k1] + c[2]*dw[k2]*g[1]*dg[k1] + 2*d[1]*dw[k2]*g[0]*g[1]*dg[k1] + 2*d[2]*dw[k2]*g[0]*g[1]*dg[k1] + 2*d[4]*dw[k2]*g[0]*g[1]*dg[k1] + d[3]*dw[k2]*std::pow(g[1],2)*dg[k1] + d[5]*dw[k2]*std::pow(g[1],2)*dg[k1] + d[6]*dw[k2]*std::pow(g[1],2)*dg[k1] + w*db[k2]*dg[k1] + 2*w*g[0]*dc[k2]*dg[k1] + w*g[1]*dc[2 + k2]*dg[k1] + w*g[1]*dc[4 + k2]*dg[k1] + 3*w*std::pow(g[0],2)*dd[k2]*dg[k1] + 2*w*g[0]*g[1]*dd[2 + k2]*dg[k1] + 2*w*g[0]*g[1]*dd[4 + k2]*dg[k1] + w*std::pow(g[1],2)*dd[6 + k2]*dg[k1] + 2*w*g[0]*g[1]*dd[8 + k2]*dg[k1] + w*std::pow(g[1],2)*dd[10 + k2]*dg[k1] + w*std::pow(g[1],2)*dd[12 + k2]*dg[k1] + dw[k2]*(g[0]*db[k1] + g[1]*db[2 + k1] + std::pow(g[0],2)*dc[k1] + g[0]*g[1]*dc[2 + k1] + g[0]*g[1]*dc[4 + k1] + std::pow(g[1],2)*dc[6 + k1] + std::pow(g[0],3)*dd[k1] + std::pow(g[0],2)*g[1]*dd[2 + k1] + std::pow(g[0],2)*g[1]*dd[4 + k1] + g[0]*std::pow(g[1],2)*dd[6 + k1] + std::pow(g[0],2)*g[1]*dd[8 + k1] + g[0]*std::pow(g[1],2)*dd[10 + k1] + g[0]*std::pow(g[1],2)*dd[12 + k1] + std::pow(g[1],3)*dd[14 + k1] + b[0]*dg[k1]) + 2*c[0]*dw[k1]*g[0]*dg[k2] + 3*d[0]*dw[k1]*std::pow(g[0],2)*dg[k2] + c[1]*dw[k1]*g[1]*dg[k2] + c[2]*dw[k1]*g[1]*dg[k2] + 2*d[1]*dw[k1]*g[0]*g[1]*dg[k2] + 2*d[2]*dw[k1]*g[0]*g[1]*dg[k2] + 2*d[4]*dw[k1]*g[0]*g[1]*dg[k2] + d[3]*dw[k1]*std::pow(g[1],2)*dg[k2] + d[5]*dw[k1]*std::pow(g[1],2)*dg[k2] + d[6]*dw[k1]*std::pow(g[1],2)*dg[k2] + w*db[k1]*dg[k2] + 2*w*g[0]*dc[k1]*dg[k2] + w*g[1]*dc[2 + k1]*dg[k2] + w*g[1]*dc[4 + k1]*dg[k2] + 3*w*std::pow(g[0],2)*dd[k1]*dg[k2] + 2*w*g[0]*g[1]*dd[2 + k1]*dg[k2] + 2*w*g[0]*g[1]*dd[4 + k1]*dg[k2] + w*std::pow(g[1],2)*dd[6 + k1]*dg[k2] + 2*w*g[0]*g[1]*dd[8 + k1]*dg[k2] + w*std::pow(g[1],2)*dd[10 + k1]*dg[k2] + w*std::pow(g[1],2)*dd[12 + k1]*dg[k2] + 2*w*c[0]*dg[k1]*dg[k2] + 6*w*d[0]*g[0]*dg[k1]*dg[k2] + 2*w*d[1]*g[1]*dg[k1]*dg[k2] + 2*w*d[2]*g[1]*dg[k1]*dg[k2] + 2*w*d[4]*g[1]*dg[k1]*dg[k2] + dw[k1]*(g[0]*db[k2] + g[1]*db[2 + k2] + std::pow(g[0],2)*dc[k2] + g[0]*g[1]*dc[2 + k2] + g[0]*g[1]*dc[4 + k2] + std::pow(g[1],2)*dc[6 + k2] + std::pow(g[0],3)*dd[k2] + std::pow(g[0],2)*g[1]*dd[2 + k2] + std::pow(g[0],2)*g[1]*dd[4 + k2] + g[0]*std::pow(g[1],2)*dd[6 + k2] + std::pow(g[0],2)*g[1]*dd[8 + k2] + g[0]*std::pow(g[1],2)*dd[10 + k2] + g[0]*std::pow(g[1],2)*dd[12 + k2] + std::pow(g[1],3)*dd[14 + k2] + b[0]*dg[k2]) + b[1]*dw[k2]*dg[2 + k1] + c[1]*dw[k2]*g[0]*dg[2 + k1] + c[2]*dw[k2]*g[0]*dg[2 + k1] + d[1]*dw[k2]*std::pow(g[0],2)*dg[2 + k1] + d[2]*dw[k2]*std::pow(g[0],2)*dg[2 + k1] + d[4]*dw[k2]*std::pow(g[0],2)*dg[2 + k1] + 2*c[3]*dw[k2]*g[1]*dg[2 + k1] + 2*d[3]*dw[k2]*g[0]*g[1]*dg[2 + k1] + 2*d[5]*dw[k2]*g[0]*g[1]*dg[2 + k1] + 2*d[6]*dw[k2]*g[0]*g[1]*dg[2 + k1] + 3*d[7]*dw[k2]*std::pow(g[1],2)*dg[2 + k1] + w*db[2 + k2]*dg[2 + k1] + w*g[0]*dc[2 + k2]*dg[2 + k1] + w*g[0]*dc[4 + k2]*dg[2 + k1] + 2*w*g[1]*dc[6 + k2]*dg[2 + k1] + w*std::pow(g[0],2)*dd[2 + k2]*dg[2 + k1] + w*std::pow(g[0],2)*dd[4 + k2]*dg[2 + k1] + 2*w*g[0]*g[1]*dd[6 + k2]*dg[2 + k1] + w*std::pow(g[0],2)*dd[8 + k2]*dg[2 + k1] + 2*w*g[0]*g[1]*dd[10 + k2]*dg[2 + k1] + 2*w*g[0]*g[1]*dd[12 + k2]*dg[2 + k1] + 3*w*std::pow(g[1],2)*dd[14 + k2]*dg[2 + k1] + w*c[1]*dg[k2]*dg[2 + k1] + w*c[2]*dg[k2]*dg[2 + k1] + 2*w*d[1]*g[0]*dg[k2]*dg[2 + k1] + 2*w*d[2]*g[0]*dg[k2]*dg[2 + k1] + 2*w*d[4]*g[0]*dg[k2]*dg[2 + k1] + 2*w*d[3]*g[1]*dg[k2]*dg[2 + k1] + 2*w*d[5]*g[1]*dg[k2]*dg[2 + k1] + 2*w*d[6]*g[1]*dg[k2]*dg[2 + k1] + (dw[k1]*(b[1] + g[0]*(c[1] + c[2] + (d[1] + d[2] + d[4])*g[0]) + 2*(c[3] + (d[3] + d[5] + d[6])*g[0])*g[1] + 3*d[7]*std::pow(g[1],2)) + w*(db[2 + k1] + g[0]*dc[2 + k1] + g[0]*dc[4 + k1] + 2*g[1]*dc[6 + k1] + std::pow(g[0],2)*(dd[2 + k1] + dd[4 + k1]) + 2*g[0]*g[1]*dd[6 + k1] + std::pow(g[0],2)*dd[8 + k1] + 2*g[0]*g[1]*dd[10 + k1] + 2*g[0]*g[1]*dd[12 + k1] + 3*std::pow(g[1],2)*dd[14 + k1] + (c[1] + c[2])*dg[k1] + 2*d[1]*g[0]*dg[k1] + 2*d[2]*g[0]*dg[k1] + 2*d[4]*g[0]*dg[k1] + 2*d[3]*g[1]*dg[k1] + 2*d[5]*g[1]*dg[k1] + 2*d[6]*g[1]*dg[k1] + 2*c[3]*dg[2 + k1] + 2*((d[3] + d[5] + d[6])*g[0] + 3*d[7]*g[1])*dg[2 + k1]))*dg[2 + k2];
      }
    }
    break;
  }
  return SymTensor(hess);
}

template<>
Dim<3>::SymTensor
evaluateRKHessian<Dim<3>>(const TableKernel<Dim<3>>& kernel,
                          const CRKOrder correctionOrder,
                          const Dim<3>::Vector& eta,
                          const Dim<3>::SymTensor& H,
                          const Dim<3>::Vector& g,
                          const Dim<3>::Tensor& dg,
                          const Dim<3>::Scalar& a,
                          const Dim<3>::Vector& b,
                          const Dim<3>::Tensor& c,
                          const Dim<3>::ThirdRankTensor& d,
                          const Dim<3>::Vector& da,
                          const Dim<3>::Tensor& db,
                          const Dim<3>::ThirdRankTensor& dc,
                          const Dim<3>::FourthRankTensor& dd,
                          const Dim<3>::Tensor& dda,
                          const Dim<3>::ThirdRankTensor& ddb,
                          const Dim<3>::FourthRankTensor& ddc,
                          const Dim<3>::FifthRankTensor& ddd) {
  typedef Dim<3>::Tensor Tensor;
  typedef Dim<3>::SymTensor SymTensor;
  typedef Dim<3>::ThirdRankTensor ThirdRankTensor;

  // Evaluate the kernel
  const auto deta = Tensor(H);
  const auto ddeta = ThirdRankTensor::zero;
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto k = kernel(eta, H);
  const auto dk = kernel.grad(eta, H);
  const auto ddk = kernel.grad2(eta, H);
  const auto w = k;
  const auto dw = Heta * dk;
  const auto ddw = Tensor((H2 - Heta2) * etaMagInv * dk + Heta2 * ddk);
  
  // Initialize values
  Tensor hess = Tensor::zero;
  const auto dim = Dim<3>::nDim;

  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::ZerothOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[3*k1 + k2] = w*dda[3*k1 + k2] + a*ddw[3*k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2];
      }
    }
    break;
  case CRKOrder::LinearOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[3*k1 + k2] = ddw[3*k1 + k2]*(a + b[0]*g[0] + b[1]*g[1] + b[2]*g[2]) + dw[k2]*(da[k1] + g[0]*db[k1] + g[1]*db[3 + k1] + g[2]*db[6 + k1] + b[0]*dg[k1] + b[1]*dg[3 + k1] + b[2]*dg[6 + k1]) + dw[k1]*(da[k2] + g[0]*db[k2] + g[1]*db[3 + k2] + g[2]*db[6 + k2] + b[0]*dg[k2] + b[1]*dg[3 + k2] + b[2]*dg[6 + k2]) + w*(dda[3*k1 + k2] + g[0]*ddb[3*k1 + k2] + g[1]*ddb[9 + 3*k1 + k2] + g[2]*ddb[18 + 3*k1 + k2] + db[k2]*dg[k1] + db[k1]*dg[k2] + db[3 + k2]*dg[3 + k1] + db[3 + k1]*dg[3 + k2] + db[6 + k2]*dg[6 + k1] + db[6 + k1]*dg[6 + k2]);
      }
    }
    break;
  case CRKOrder::QuadraticOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[3*k1 + k2] = w*dda[3*k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2] + ddw[3*k1 + k2]*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + g[1]*(b[1] + (c[1] + c[3])*g[0] + c[4]*g[1]) + (b[2] + (c[2] + c[6])*g[0] + (c[5] + c[7])*g[1])*g[2] + c[8]*std::pow(g[2],2)) + w*g[0]*ddb[3*k1 + k2] + w*g[1]*ddb[9 + 3*k1 + k2] + w*g[2]*ddb[18 + 3*k1 + k2] + w*std::pow(g[0],2)*ddc[3*k1 + k2] + w*g[0]*g[1]*ddc[9 + 3*k1 + k2] + w*g[0]*g[2]*ddc[18 + 3*k1 + k2] + w*g[0]*g[1]*ddc[27 + 3*k1 + k2] + w*std::pow(g[1],2)*ddc[36 + 3*k1 + k2] + w*g[1]*g[2]*ddc[45 + 3*k1 + k2] + w*g[0]*g[2]*ddc[54 + 3*k1 + k2] + w*g[1]*g[2]*ddc[63 + 3*k1 + k2] + w*std::pow(g[2],2)*ddc[72 + 3*k1 + k2] + 2*c[0]*dw[k2]*g[0]*dg[k1] + c[1]*dw[k2]*g[1]*dg[k1] + c[3]*dw[k2]*g[1]*dg[k1] + c[2]*dw[k2]*g[2]*dg[k1] + c[6]*dw[k2]*g[2]*dg[k1] + w*db[k2]*dg[k1] + 2*w*g[0]*dc[k2]*dg[k1] + w*g[1]*dc[3 + k2]*dg[k1] + w*g[2]*dc[6 + k2]*dg[k1] + w*g[1]*dc[9 + k2]*dg[k1] + w*g[2]*dc[18 + k2]*dg[k1] + dw[k2]*(g[0]*db[k1] + g[1]*db[3 + k1] + g[2]*db[6 + k1] + std::pow(g[0],2)*dc[k1] + g[0]*g[1]*dc[3 + k1] + g[0]*g[2]*dc[6 + k1] + g[0]*g[1]*dc[9 + k1] + std::pow(g[1],2)*dc[12 + k1] + g[1]*g[2]*dc[15 + k1] + g[0]*g[2]*dc[18 + k1] + g[1]*g[2]*dc[21 + k1] + std::pow(g[2],2)*dc[24 + k1] + b[0]*dg[k1]) + 2*c[0]*dw[k1]*g[0]*dg[k2] + c[1]*dw[k1]*g[1]*dg[k2] + c[3]*dw[k1]*g[1]*dg[k2] + c[2]*dw[k1]*g[2]*dg[k2] + c[6]*dw[k1]*g[2]*dg[k2] + w*db[k1]*dg[k2] + 2*w*g[0]*dc[k1]*dg[k2] + w*g[1]*dc[3 + k1]*dg[k2] + w*g[2]*dc[6 + k1]*dg[k2] + w*g[1]*dc[9 + k1]*dg[k2] + w*g[2]*dc[18 + k1]*dg[k2] + 2*w*c[0]*dg[k1]*dg[k2] + dw[k1]*(g[0]*db[k2] + g[1]*db[3 + k2] + g[2]*db[6 + k2] + std::pow(g[0],2)*dc[k2] + g[0]*g[1]*dc[3 + k2] + g[0]*g[2]*dc[6 + k2] + g[0]*g[1]*dc[9 + k2] + std::pow(g[1],2)*dc[12 + k2] + g[1]*g[2]*dc[15 + k2] + g[0]*g[2]*dc[18 + k2] + g[1]*g[2]*dc[21 + k2] + std::pow(g[2],2)*dc[24 + k2] + b[0]*dg[k2]) + b[1]*dw[k2]*dg[3 + k1] + c[1]*dw[k2]*g[0]*dg[3 + k1] + c[3]*dw[k2]*g[0]*dg[3 + k1] + 2*c[4]*dw[k2]*g[1]*dg[3 + k1] + c[5]*dw[k2]*g[2]*dg[3 + k1] + c[7]*dw[k2]*g[2]*dg[3 + k1] + w*db[3 + k2]*dg[3 + k1] + w*g[0]*dc[3 + k2]*dg[3 + k1] + w*g[0]*dc[9 + k2]*dg[3 + k1] + 2*w*g[1]*dc[12 + k2]*dg[3 + k1] + w*g[2]*dc[15 + k2]*dg[3 + k1] + w*g[2]*dc[21 + k2]*dg[3 + k1] + w*c[1]*dg[k2]*dg[3 + k1] + w*c[3]*dg[k2]*dg[3 + k1] + b[1]*dw[k1]*dg[3 + k2] + c[1]*dw[k1]*g[0]*dg[3 + k2] + c[3]*dw[k1]*g[0]*dg[3 + k2] + 2*c[4]*dw[k1]*g[1]*dg[3 + k2] + c[5]*dw[k1]*g[2]*dg[3 + k2] + c[7]*dw[k1]*g[2]*dg[3 + k2] + w*db[3 + k1]*dg[3 + k2] + w*g[0]*dc[3 + k1]*dg[3 + k2] + w*g[0]*dc[9 + k1]*dg[3 + k2] + 2*w*g[1]*dc[12 + k1]*dg[3 + k2] + w*g[2]*dc[15 + k1]*dg[3 + k2] + w*g[2]*dc[21 + k1]*dg[3 + k2] + w*c[1]*dg[k1]*dg[3 + k2] + w*c[3]*dg[k1]*dg[3 + k2] + 2*w*c[4]*dg[3 + k1]*dg[3 + k2] + b[2]*dw[k2]*dg[6 + k1] + c[2]*dw[k2]*g[0]*dg[6 + k1] + c[6]*dw[k2]*g[0]*dg[6 + k1] + c[5]*dw[k2]*g[1]*dg[6 + k1] + c[7]*dw[k2]*g[1]*dg[6 + k1] + 2*c[8]*dw[k2]*g[2]*dg[6 + k1] + w*db[6 + k2]*dg[6 + k1] + w*g[0]*dc[6 + k2]*dg[6 + k1] + w*g[1]*dc[15 + k2]*dg[6 + k1] + w*g[0]*dc[18 + k2]*dg[6 + k1] + w*g[1]*dc[21 + k2]*dg[6 + k1] + 2*w*g[2]*dc[24 + k2]*dg[6 + k1] + w*c[2]*dg[k2]*dg[6 + k1] + w*c[6]*dg[k2]*dg[6 + k1] + w*c[5]*dg[3 + k2]*dg[6 + k1] + w*c[7]*dg[3 + k2]*dg[6 + k1] + (dw[k1]*(b[2] + (c[2] + c[6])*g[0] + (c[5] + c[7])*g[1] + 2*c[8]*g[2]) + w*(db[6 + k1] + g[0]*(dc[6 + k1] + dc[18 + k1]) + g[1]*(dc[15 + k1] + dc[21 + k1]) + 2*g[2]*dc[24 + k1] + (c[2] + c[6])*dg[k1] + (c[5] + c[7])*dg[3 + k1] + 2*c[8]*dg[6 + k1]))*dg[6 + k2];
      }
    }
    break;
  case CRKOrder::CubicOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[3*k1 + k2] = w*dda[3*k1 + k2] + da[k2]*dw[k1] + da[k1]*dw[k2] + ddw[3*k1 + k2]*(a + b[0]*g[0] + c[0]*std::pow(g[0],2) + d[0]*std::pow(g[0],3) + b[1]*g[1] + c[1]*g[0]*g[1] + c[3]*g[0]*g[1] + d[1]*std::pow(g[0],2)*g[1] + d[3]*std::pow(g[0],2)*g[1] + d[9]*std::pow(g[0],2)*g[1] + c[4]*std::pow(g[1],2) + d[4]*g[0]*std::pow(g[1],2) + d[10]*g[0]*std::pow(g[1],2) + d[12]*g[0]*std::pow(g[1],2) + d[13]*std::pow(g[1],3) + (b[2] + g[0]*(c[2] + c[6] + (d[2] + d[6] + d[18])*g[0]) + (c[5] + c[7] + (d[5] + d[7] + d[11] + d[15] + d[19] + d[21])*g[0])*g[1] + (d[14] + d[16] + d[22])*std::pow(g[1],2))*g[2] + (c[8] + (d[8] + d[20] + d[24])*g[0] + (d[17] + d[23] + d[25])*g[1])*std::pow(g[2],2) + d[26]*std::pow(g[2],3)) + w*g[0]*ddb[3*k1 + k2] + w*g[1]*ddb[9 + 3*k1 + k2] + w*g[2]*ddb[18 + 3*k1 + k2] + w*std::pow(g[0],2)*ddc[3*k1 + k2] + w*g[0]*g[1]*ddc[9 + 3*k1 + k2] + w*g[0]*g[2]*ddc[18 + 3*k1 + k2] + w*g[0]*g[1]*ddc[27 + 3*k1 + k2] + w*std::pow(g[1],2)*ddc[36 + 3*k1 + k2] + w*g[1]*g[2]*ddc[45 + 3*k1 + k2] + w*g[0]*g[2]*ddc[54 + 3*k1 + k2] + w*g[1]*g[2]*ddc[63 + 3*k1 + k2] + w*std::pow(g[2],2)*ddc[72 + 3*k1 + k2] + w*std::pow(g[0],3)*ddd[3*k1 + k2] + w*std::pow(g[0],2)*g[1]*ddd[9 + 3*k1 + k2] + w*std::pow(g[0],2)*g[2]*ddd[18 + 3*k1 + k2] + w*std::pow(g[0],2)*g[1]*ddd[27 + 3*k1 + k2] + w*g[0]*std::pow(g[1],2)*ddd[36 + 3*k1 + k2] + w*g[0]*g[1]*g[2]*ddd[45 + 3*k1 + k2] + w*std::pow(g[0],2)*g[2]*ddd[54 + 3*k1 + k2] + w*g[0]*g[1]*g[2]*ddd[63 + 3*k1 + k2] + w*g[0]*std::pow(g[2],2)*ddd[72 + 3*k1 + k2] + w*std::pow(g[0],2)*g[1]*ddd[81 + 3*k1 + k2] + w*g[0]*std::pow(g[1],2)*ddd[90 + 3*k1 + k2] + w*g[0]*g[1]*g[2]*ddd[99 + 3*k1 + k2] + w*g[0]*std::pow(g[1],2)*ddd[108 + 3*k1 + k2] + w*std::pow(g[1],3)*ddd[117 + 3*k1 + k2] + w*std::pow(g[1],2)*g[2]*ddd[126 + 3*k1 + k2] + w*g[0]*g[1]*g[2]*ddd[135 + 3*k1 + k2] + w*std::pow(g[1],2)*g[2]*ddd[144 + 3*k1 + k2] + w*g[1]*std::pow(g[2],2)*ddd[153 + 3*k1 + k2] + w*std::pow(g[0],2)*g[2]*ddd[162 + 3*k1 + k2] + w*g[0]*g[1]*g[2]*ddd[171 + 3*k1 + k2] + w*g[0]*std::pow(g[2],2)*ddd[180 + 3*k1 + k2] + w*g[0]*g[1]*g[2]*ddd[189 + 3*k1 + k2] + w*std::pow(g[1],2)*g[2]*ddd[198 + 3*k1 + k2] + w*g[1]*std::pow(g[2],2)*ddd[207 + 3*k1 + k2] + w*g[0]*std::pow(g[2],2)*ddd[216 + 3*k1 + k2] + w*g[1]*std::pow(g[2],2)*ddd[225 + 3*k1 + k2] + w*std::pow(g[2],3)*ddd[234 + 3*k1 + k2] + 2*c[0]*dw[k2]*g[0]*dg[k1] + 3*d[0]*dw[k2]*std::pow(g[0],2)*dg[k1] + c[1]*dw[k2]*g[1]*dg[k1] + c[3]*dw[k2]*g[1]*dg[k1] + 2*d[1]*dw[k2]*g[0]*g[1]*dg[k1] + 2*d[3]*dw[k2]*g[0]*g[1]*dg[k1] + 2*d[9]*dw[k2]*g[0]*g[1]*dg[k1] + d[4]*dw[k2]*std::pow(g[1],2)*dg[k1] + d[10]*dw[k2]*std::pow(g[1],2)*dg[k1] + d[12]*dw[k2]*std::pow(g[1],2)*dg[k1] + c[2]*dw[k2]*g[2]*dg[k1] + c[6]*dw[k2]*g[2]*dg[k1] + 2*d[2]*dw[k2]*g[0]*g[2]*dg[k1] + 2*d[6]*dw[k2]*g[0]*g[2]*dg[k1] + 2*d[18]*dw[k2]*g[0]*g[2]*dg[k1] + d[5]*dw[k2]*g[1]*g[2]*dg[k1] + d[7]*dw[k2]*g[1]*g[2]*dg[k1] + d[11]*dw[k2]*g[1]*g[2]*dg[k1] + d[15]*dw[k2]*g[1]*g[2]*dg[k1] + d[19]*dw[k2]*g[1]*g[2]*dg[k1] + d[21]*dw[k2]*g[1]*g[2]*dg[k1] + d[8]*dw[k2]*std::pow(g[2],2)*dg[k1] + d[20]*dw[k2]*std::pow(g[2],2)*dg[k1] + d[24]*dw[k2]*std::pow(g[2],2)*dg[k1] + w*db[k2]*dg[k1] + 2*w*g[0]*dc[k2]*dg[k1] + w*g[1]*dc[3 + k2]*dg[k1] + w*g[2]*dc[6 + k2]*dg[k1] + w*g[1]*dc[9 + k2]*dg[k1] + w*g[2]*dc[18 + k2]*dg[k1] + 3*w*std::pow(g[0],2)*dd[k2]*dg[k1] + 2*w*g[0]*g[1]*dd[3 + k2]*dg[k1] + 2*w*g[0]*g[2]*dd[6 + k2]*dg[k1] + 2*w*g[0]*g[1]*dd[9 + k2]*dg[k1] + w*std::pow(g[1],2)*dd[12 + k2]*dg[k1] + w*g[1]*g[2]*dd[15 + k2]*dg[k1] + 2*w*g[0]*g[2]*dd[18 + k2]*dg[k1] + w*g[1]*g[2]*dd[21 + k2]*dg[k1] + w*std::pow(g[2],2)*dd[24 + k2]*dg[k1] + 2*w*g[0]*g[1]*dd[27 + k2]*dg[k1] + w*std::pow(g[1],2)*dd[30 + k2]*dg[k1] + w*g[1]*g[2]*dd[33 + k2]*dg[k1] + w*std::pow(g[1],2)*dd[36 + k2]*dg[k1] + w*g[1]*g[2]*dd[45 + k2]*dg[k1] + 2*w*g[0]*g[2]*dd[54 + k2]*dg[k1] + w*g[1]*g[2]*dd[57 + k2]*dg[k1] + w*std::pow(g[2],2)*dd[60 + k2]*dg[k1] + w*g[1]*g[2]*dd[63 + k2]*dg[k1] + w*std::pow(g[2],2)*dd[72 + k2]*dg[k1] + dw[k2]*(g[0]*db[k1] + g[1]*db[3 + k1] + g[2]*db[6 + k1] + std::pow(g[0],2)*dc[k1] + g[0]*g[1]*dc[3 + k1] + g[0]*g[2]*dc[6 + k1] + g[0]*g[1]*dc[9 + k1] + std::pow(g[1],2)*dc[12 + k1] + g[1]*g[2]*dc[15 + k1] + g[0]*g[2]*dc[18 + k1] + g[1]*g[2]*dc[21 + k1] + std::pow(g[2],2)*dc[24 + k1] + std::pow(g[0],3)*dd[k1] + std::pow(g[0],2)*g[1]*dd[3 + k1] + std::pow(g[0],2)*g[2]*dd[6 + k1] + std::pow(g[0],2)*g[1]*dd[9 + k1] + g[0]*std::pow(g[1],2)*dd[12 + k1] + g[0]*g[1]*g[2]*dd[15 + k1] + std::pow(g[0],2)*g[2]*dd[18 + k1] + g[0]*g[1]*g[2]*dd[21 + k1] + g[0]*std::pow(g[2],2)*dd[24 + k1] + std::pow(g[0],2)*g[1]*dd[27 + k1] + g[0]*std::pow(g[1],2)*dd[30 + k1] + g[0]*g[1]*g[2]*dd[33 + k1] + g[0]*std::pow(g[1],2)*dd[36 + k1] + std::pow(g[1],3)*dd[39 + k1] + std::pow(g[1],2)*g[2]*dd[42 + k1] + g[0]*g[1]*g[2]*dd[45 + k1] + std::pow(g[1],2)*g[2]*dd[48 + k1] + g[1]*std::pow(g[2],2)*dd[51 + k1] + std::pow(g[0],2)*g[2]*dd[54 + k1] + g[0]*g[1]*g[2]*dd[57 + k1] + g[0]*std::pow(g[2],2)*dd[60 + k1] + g[0]*g[1]*g[2]*dd[63 + k1] + std::pow(g[1],2)*g[2]*dd[66 + k1] + g[1]*std::pow(g[2],2)*dd[69 + k1] + g[0]*std::pow(g[2],2)*dd[72 + k1] + g[1]*std::pow(g[2],2)*dd[75 + k1] + std::pow(g[2],3)*dd[78 + k1] + b[0]*dg[k1]) + 2*c[0]*dw[k1]*g[0]*dg[k2] + 3*d[0]*dw[k1]*std::pow(g[0],2)*dg[k2] + c[1]*dw[k1]*g[1]*dg[k2] + c[3]*dw[k1]*g[1]*dg[k2] + 2*d[1]*dw[k1]*g[0]*g[1]*dg[k2] + 2*d[3]*dw[k1]*g[0]*g[1]*dg[k2] + 2*d[9]*dw[k1]*g[0]*g[1]*dg[k2] + d[4]*dw[k1]*std::pow(g[1],2)*dg[k2] + d[10]*dw[k1]*std::pow(g[1],2)*dg[k2] + d[12]*dw[k1]*std::pow(g[1],2)*dg[k2] + c[2]*dw[k1]*g[2]*dg[k2] + c[6]*dw[k1]*g[2]*dg[k2] + 2*d[2]*dw[k1]*g[0]*g[2]*dg[k2] + 2*d[6]*dw[k1]*g[0]*g[2]*dg[k2] + 2*d[18]*dw[k1]*g[0]*g[2]*dg[k2] + d[5]*dw[k1]*g[1]*g[2]*dg[k2] + d[7]*dw[k1]*g[1]*g[2]*dg[k2] + d[11]*dw[k1]*g[1]*g[2]*dg[k2] + d[15]*dw[k1]*g[1]*g[2]*dg[k2] + d[19]*dw[k1]*g[1]*g[2]*dg[k2] + d[21]*dw[k1]*g[1]*g[2]*dg[k2] + d[8]*dw[k1]*std::pow(g[2],2)*dg[k2] + d[20]*dw[k1]*std::pow(g[2],2)*dg[k2] + d[24]*dw[k1]*std::pow(g[2],2)*dg[k2] + w*db[k1]*dg[k2] + 2*w*g[0]*dc[k1]*dg[k2] + w*g[1]*dc[3 + k1]*dg[k2] + w*g[2]*dc[6 + k1]*dg[k2] + w*g[1]*dc[9 + k1]*dg[k2] + w*g[2]*dc[18 + k1]*dg[k2] + 3*w*std::pow(g[0],2)*dd[k1]*dg[k2] + 2*w*g[0]*g[1]*dd[3 + k1]*dg[k2] + 2*w*g[0]*g[2]*dd[6 + k1]*dg[k2] + 2*w*g[0]*g[1]*dd[9 + k1]*dg[k2] + w*std::pow(g[1],2)*dd[12 + k1]*dg[k2] + w*g[1]*g[2]*dd[15 + k1]*dg[k2] + 2*w*g[0]*g[2]*dd[18 + k1]*dg[k2] + w*g[1]*g[2]*dd[21 + k1]*dg[k2] + w*std::pow(g[2],2)*dd[24 + k1]*dg[k2] + 2*w*g[0]*g[1]*dd[27 + k1]*dg[k2] + w*std::pow(g[1],2)*dd[30 + k1]*dg[k2] + w*g[1]*g[2]*dd[33 + k1]*dg[k2] + w*std::pow(g[1],2)*dd[36 + k1]*dg[k2] + w*g[1]*g[2]*dd[45 + k1]*dg[k2] + 2*w*g[0]*g[2]*dd[54 + k1]*dg[k2] + w*g[1]*g[2]*dd[57 + k1]*dg[k2] + w*std::pow(g[2],2)*dd[60 + k1]*dg[k2] + w*g[1]*g[2]*dd[63 + k1]*dg[k2] + w*std::pow(g[2],2)*dd[72 + k1]*dg[k2] + 2*w*c[0]*dg[k1]*dg[k2] + 6*w*d[0]*g[0]*dg[k1]*dg[k2] + 2*w*d[1]*g[1]*dg[k1]*dg[k2] + 2*w*d[3]*g[1]*dg[k1]*dg[k2] + 2*w*d[9]*g[1]*dg[k1]*dg[k2] + 2*w*d[2]*g[2]*dg[k1]*dg[k2] + 2*w*d[6]*g[2]*dg[k1]*dg[k2] + 2*w*d[18]*g[2]*dg[k1]*dg[k2] + dw[k1]*(g[0]*db[k2] + g[1]*db[3 + k2] + g[2]*db[6 + k2] + std::pow(g[0],2)*dc[k2] + g[0]*g[1]*dc[3 + k2] + g[0]*g[2]*dc[6 + k2] + g[0]*g[1]*dc[9 + k2] + std::pow(g[1],2)*dc[12 + k2] + g[1]*g[2]*dc[15 + k2] + g[0]*g[2]*dc[18 + k2] + g[1]*g[2]*dc[21 + k2] + std::pow(g[2],2)*dc[24 + k2] + std::pow(g[0],3)*dd[k2] + std::pow(g[0],2)*g[1]*dd[3 + k2] + std::pow(g[0],2)*g[2]*dd[6 + k2] + std::pow(g[0],2)*g[1]*dd[9 + k2] + g[0]*std::pow(g[1],2)*dd[12 + k2] + g[0]*g[1]*g[2]*dd[15 + k2] + std::pow(g[0],2)*g[2]*dd[18 + k2] + g[0]*g[1]*g[2]*dd[21 + k2] + g[0]*std::pow(g[2],2)*dd[24 + k2] + std::pow(g[0],2)*g[1]*dd[27 + k2] + g[0]*std::pow(g[1],2)*dd[30 + k2] + g[0]*g[1]*g[2]*dd[33 + k2] + g[0]*std::pow(g[1],2)*dd[36 + k2] + std::pow(g[1],3)*dd[39 + k2] + std::pow(g[1],2)*g[2]*dd[42 + k2] + g[0]*g[1]*g[2]*dd[45 + k2] + std::pow(g[1],2)*g[2]*dd[48 + k2] + g[1]*std::pow(g[2],2)*dd[51 + k2] + std::pow(g[0],2)*g[2]*dd[54 + k2] + g[0]*g[1]*g[2]*dd[57 + k2] + g[0]*std::pow(g[2],2)*dd[60 + k2] + g[0]*g[1]*g[2]*dd[63 + k2] + std::pow(g[1],2)*g[2]*dd[66 + k2] + g[1]*std::pow(g[2],2)*dd[69 + k2] + g[0]*std::pow(g[2],2)*dd[72 + k2] + g[1]*std::pow(g[2],2)*dd[75 + k2] + std::pow(g[2],3)*dd[78 + k2] + b[0]*dg[k2]) + b[1]*dw[k2]*dg[3 + k1] + c[1]*dw[k2]*g[0]*dg[3 + k1] + c[3]*dw[k2]*g[0]*dg[3 + k1] + d[1]*dw[k2]*std::pow(g[0],2)*dg[3 + k1] + d[3]*dw[k2]*std::pow(g[0],2)*dg[3 + k1] + d[9]*dw[k2]*std::pow(g[0],2)*dg[3 + k1] + 2*c[4]*dw[k2]*g[1]*dg[3 + k1] + 2*d[4]*dw[k2]*g[0]*g[1]*dg[3 + k1] + 2*d[10]*dw[k2]*g[0]*g[1]*dg[3 + k1] + 2*d[12]*dw[k2]*g[0]*g[1]*dg[3 + k1] + 3*d[13]*dw[k2]*std::pow(g[1],2)*dg[3 + k1] + c[5]*dw[k2]*g[2]*dg[3 + k1] + c[7]*dw[k2]*g[2]*dg[3 + k1] + d[5]*dw[k2]*g[0]*g[2]*dg[3 + k1] + d[7]*dw[k2]*g[0]*g[2]*dg[3 + k1] + d[11]*dw[k2]*g[0]*g[2]*dg[3 + k1] + d[15]*dw[k2]*g[0]*g[2]*dg[3 + k1] + d[19]*dw[k2]*g[0]*g[2]*dg[3 + k1] + d[21]*dw[k2]*g[0]*g[2]*dg[3 + k1] + 2*d[14]*dw[k2]*g[1]*g[2]*dg[3 + k1] + 2*d[16]*dw[k2]*g[1]*g[2]*dg[3 + k1] + 2*d[22]*dw[k2]*g[1]*g[2]*dg[3 + k1] + d[17]*dw[k2]*std::pow(g[2],2)*dg[3 + k1] + d[23]*dw[k2]*std::pow(g[2],2)*dg[3 + k1] + d[25]*dw[k2]*std::pow(g[2],2)*dg[3 + k1] + w*db[3 + k2]*dg[3 + k1] + w*g[0]*dc[3 + k2]*dg[3 + k1] + w*g[0]*dc[9 + k2]*dg[3 + k1] + 2*w*g[1]*dc[12 + k2]*dg[3 + k1] + w*g[2]*dc[15 + k2]*dg[3 + k1] + w*g[2]*dc[21 + k2]*dg[3 + k1] + w*std::pow(g[0],2)*dd[3 + k2]*dg[3 + k1] + w*std::pow(g[0],2)*dd[9 + k2]*dg[3 + k1] + 2*w*g[0]*g[1]*dd[12 + k2]*dg[3 + k1] + w*g[0]*g[2]*dd[15 + k2]*dg[3 + k1] + w*g[0]*g[2]*dd[21 + k2]*dg[3 + k1] + w*std::pow(g[0],2)*dd[27 + k2]*dg[3 + k1] + 2*w*g[0]*g[1]*dd[30 + k2]*dg[3 + k1] + w*g[0]*g[2]*dd[33 + k2]*dg[3 + k1] + 2*w*g[0]*g[1]*dd[36 + k2]*dg[3 + k1] + 3*w*std::pow(g[1],2)*dd[39 + k2]*dg[3 + k1] + 2*w*g[1]*g[2]*dd[42 + k2]*dg[3 + k1] + w*g[0]*g[2]*dd[45 + k2]*dg[3 + k1] + 2*w*g[1]*g[2]*dd[48 + k2]*dg[3 + k1] + w*std::pow(g[2],2)*dd[51 + k2]*dg[3 + k1] + w*g[0]*g[2]*dd[57 + k2]*dg[3 + k1] + w*g[0]*g[2]*dd[63 + k2]*dg[3 + k1] + 2*w*g[1]*g[2]*dd[66 + k2]*dg[3 + k1] + w*std::pow(g[2],2)*dd[69 + k2]*dg[3 + k1] + w*std::pow(g[2],2)*dd[75 + k2]*dg[3 + k1] + w*c[1]*dg[k2]*dg[3 + k1] + w*c[3]*dg[k2]*dg[3 + k1] + 2*w*d[1]*g[0]*dg[k2]*dg[3 + k1] + 2*w*d[3]*g[0]*dg[k2]*dg[3 + k1] + 2*w*d[9]*g[0]*dg[k2]*dg[3 + k1] + 2*w*d[4]*g[1]*dg[k2]*dg[3 + k1] + 2*w*d[10]*g[1]*dg[k2]*dg[3 + k1] + 2*w*d[12]*g[1]*dg[k2]*dg[3 + k1] + w*d[5]*g[2]*dg[k2]*dg[3 + k1] + w*d[7]*g[2]*dg[k2]*dg[3 + k1] + w*d[11]*g[2]*dg[k2]*dg[3 + k1] + w*d[15]*g[2]*dg[k2]*dg[3 + k1] + w*d[19]*g[2]*dg[k2]*dg[3 + k1] + w*d[21]*g[2]*dg[k2]*dg[3 + k1] + b[1]*dw[k1]*dg[3 + k2] + c[1]*dw[k1]*g[0]*dg[3 + k2] + c[3]*dw[k1]*g[0]*dg[3 + k2] + d[1]*dw[k1]*std::pow(g[0],2)*dg[3 + k2] + d[3]*dw[k1]*std::pow(g[0],2)*dg[3 + k2] + d[9]*dw[k1]*std::pow(g[0],2)*dg[3 + k2] + 2*c[4]*dw[k1]*g[1]*dg[3 + k2] + 2*d[4]*dw[k1]*g[0]*g[1]*dg[3 + k2] + 2*d[10]*dw[k1]*g[0]*g[1]*dg[3 + k2] + 2*d[12]*dw[k1]*g[0]*g[1]*dg[3 + k2] + 3*d[13]*dw[k1]*std::pow(g[1],2)*dg[3 + k2] + c[5]*dw[k1]*g[2]*dg[3 + k2] + c[7]*dw[k1]*g[2]*dg[3 + k2] + d[5]*dw[k1]*g[0]*g[2]*dg[3 + k2] + d[7]*dw[k1]*g[0]*g[2]*dg[3 + k2] + d[11]*dw[k1]*g[0]*g[2]*dg[3 + k2] + d[15]*dw[k1]*g[0]*g[2]*dg[3 + k2] + d[19]*dw[k1]*g[0]*g[2]*dg[3 + k2] + d[21]*dw[k1]*g[0]*g[2]*dg[3 + k2] + 2*d[14]*dw[k1]*g[1]*g[2]*dg[3 + k2] + 2*d[16]*dw[k1]*g[1]*g[2]*dg[3 + k2] + 2*d[22]*dw[k1]*g[1]*g[2]*dg[3 + k2] + d[17]*dw[k1]*std::pow(g[2],2)*dg[3 + k2] + d[23]*dw[k1]*std::pow(g[2],2)*dg[3 + k2] + d[25]*dw[k1]*std::pow(g[2],2)*dg[3 + k2] + w*db[3 + k1]*dg[3 + k2] + w*g[0]*dc[3 + k1]*dg[3 + k2] + w*g[0]*dc[9 + k1]*dg[3 + k2] + 2*w*g[1]*dc[12 + k1]*dg[3 + k2] + w*g[2]*dc[15 + k1]*dg[3 + k2] + w*g[2]*dc[21 + k1]*dg[3 + k2] + w*std::pow(g[0],2)*dd[3 + k1]*dg[3 + k2] + w*std::pow(g[0],2)*dd[9 + k1]*dg[3 + k2] + 2*w*g[0]*g[1]*dd[12 + k1]*dg[3 + k2] + w*g[0]*g[2]*dd[15 + k1]*dg[3 + k2] + w*g[0]*g[2]*dd[21 + k1]*dg[3 + k2] + w*std::pow(g[0],2)*dd[27 + k1]*dg[3 + k2] + 2*w*g[0]*g[1]*dd[30 + k1]*dg[3 + k2] + w*g[0]*g[2]*dd[33 + k1]*dg[3 + k2] + 2*w*g[0]*g[1]*dd[36 + k1]*dg[3 + k2] + 3*w*std::pow(g[1],2)*dd[39 + k1]*dg[3 + k2] + 2*w*g[1]*g[2]*dd[42 + k1]*dg[3 + k2] + w*g[0]*g[2]*dd[45 + k1]*dg[3 + k2] + 2*w*g[1]*g[2]*dd[48 + k1]*dg[3 + k2] + w*std::pow(g[2],2)*dd[51 + k1]*dg[3 + k2] + w*g[0]*g[2]*dd[57 + k1]*dg[3 + k2] + w*g[0]*g[2]*dd[63 + k1]*dg[3 + k2] + 2*w*g[1]*g[2]*dd[66 + k1]*dg[3 + k2] + w*std::pow(g[2],2)*dd[69 + k1]*dg[3 + k2] + w*std::pow(g[2],2)*dd[75 + k1]*dg[3 + k2] + w*c[1]*dg[k1]*dg[3 + k2] + w*c[3]*dg[k1]*dg[3 + k2] + 2*w*d[1]*g[0]*dg[k1]*dg[3 + k2] + 2*w*d[3]*g[0]*dg[k1]*dg[3 + k2] + 2*w*d[9]*g[0]*dg[k1]*dg[3 + k2] + 2*w*d[4]*g[1]*dg[k1]*dg[3 + k2] + 2*w*d[10]*g[1]*dg[k1]*dg[3 + k2] + 2*w*d[12]*g[1]*dg[k1]*dg[3 + k2] + w*d[5]*g[2]*dg[k1]*dg[3 + k2] + w*d[7]*g[2]*dg[k1]*dg[3 + k2] + w*d[11]*g[2]*dg[k1]*dg[3 + k2] + w*d[15]*g[2]*dg[k1]*dg[3 + k2] + w*d[19]*g[2]*dg[k1]*dg[3 + k2] + w*d[21]*g[2]*dg[k1]*dg[3 + k2] + 2*w*c[4]*dg[3 + k1]*dg[3 + k2] + 2*w*d[4]*g[0]*dg[3 + k1]*dg[3 + k2] + 2*w*d[10]*g[0]*dg[3 + k1]*dg[3 + k2] + 2*w*d[12]*g[0]*dg[3 + k1]*dg[3 + k2] + 6*w*d[13]*g[1]*dg[3 + k1]*dg[3 + k2] + 2*w*d[14]*g[2]*dg[3 + k1]*dg[3 + k2] + 2*w*d[16]*g[2]*dg[3 + k1]*dg[3 + k2] + 2*w*d[22]*g[2]*dg[3 + k1]*dg[3 + k2] + b[2]*dw[k2]*dg[6 + k1] + c[2]*dw[k2]*g[0]*dg[6 + k1] + c[6]*dw[k2]*g[0]*dg[6 + k1] + d[2]*dw[k2]*std::pow(g[0],2)*dg[6 + k1] + d[6]*dw[k2]*std::pow(g[0],2)*dg[6 + k1] + d[18]*dw[k2]*std::pow(g[0],2)*dg[6 + k1] + c[5]*dw[k2]*g[1]*dg[6 + k1] + c[7]*dw[k2]*g[1]*dg[6 + k1] + d[5]*dw[k2]*g[0]*g[1]*dg[6 + k1] + d[7]*dw[k2]*g[0]*g[1]*dg[6 + k1] + d[11]*dw[k2]*g[0]*g[1]*dg[6 + k1] + d[15]*dw[k2]*g[0]*g[1]*dg[6 + k1] + d[19]*dw[k2]*g[0]*g[1]*dg[6 + k1] + d[21]*dw[k2]*g[0]*g[1]*dg[6 + k1] + d[14]*dw[k2]*std::pow(g[1],2)*dg[6 + k1] + d[16]*dw[k2]*std::pow(g[1],2)*dg[6 + k1] + d[22]*dw[k2]*std::pow(g[1],2)*dg[6 + k1] + 2*c[8]*dw[k2]*g[2]*dg[6 + k1] + 2*d[8]*dw[k2]*g[0]*g[2]*dg[6 + k1] + 2*d[20]*dw[k2]*g[0]*g[2]*dg[6 + k1] + 2*d[24]*dw[k2]*g[0]*g[2]*dg[6 + k1] + 2*d[17]*dw[k2]*g[1]*g[2]*dg[6 + k1] + 2*d[23]*dw[k2]*g[1]*g[2]*dg[6 + k1] + 2*d[25]*dw[k2]*g[1]*g[2]*dg[6 + k1] + 3*d[26]*dw[k2]*std::pow(g[2],2)*dg[6 + k1] + w*db[6 + k2]*dg[6 + k1] + w*g[0]*dc[6 + k2]*dg[6 + k1] + w*g[1]*dc[15 + k2]*dg[6 + k1] + w*g[0]*dc[18 + k2]*dg[6 + k1] + w*g[1]*dc[21 + k2]*dg[6 + k1] + 2*w*g[2]*dc[24 + k2]*dg[6 + k1] + w*std::pow(g[0],2)*dd[6 + k2]*dg[6 + k1] + w*g[0]*g[1]*dd[15 + k2]*dg[6 + k1] + w*std::pow(g[0],2)*dd[18 + k2]*dg[6 + k1] + w*g[0]*g[1]*dd[21 + k2]*dg[6 + k1] + 2*w*g[0]*g[2]*dd[24 + k2]*dg[6 + k1] + w*g[0]*g[1]*dd[33 + k2]*dg[6 + k1] + w*std::pow(g[1],2)*dd[42 + k2]*dg[6 + k1] + w*g[0]*g[1]*dd[45 + k2]*dg[6 + k1] + w*std::pow(g[1],2)*dd[48 + k2]*dg[6 + k1] + 2*w*g[1]*g[2]*dd[51 + k2]*dg[6 + k1] + w*std::pow(g[0],2)*dd[54 + k2]*dg[6 + k1] + w*g[0]*g[1]*dd[57 + k2]*dg[6 + k1] + 2*w*g[0]*g[2]*dd[60 + k2]*dg[6 + k1] + w*g[0]*g[1]*dd[63 + k2]*dg[6 + k1] + w*std::pow(g[1],2)*dd[66 + k2]*dg[6 + k1] + 2*w*g[1]*g[2]*dd[69 + k2]*dg[6 + k1] + 2*w*g[0]*g[2]*dd[72 + k2]*dg[6 + k1] + 2*w*g[1]*g[2]*dd[75 + k2]*dg[6 + k1] + 3*w*std::pow(g[2],2)*dd[78 + k2]*dg[6 + k1] + w*c[2]*dg[k2]*dg[6 + k1] + w*c[6]*dg[k2]*dg[6 + k1] + 2*w*d[2]*g[0]*dg[k2]*dg[6 + k1] + 2*w*d[6]*g[0]*dg[k2]*dg[6 + k1] + 2*w*d[18]*g[0]*dg[k2]*dg[6 + k1] + w*d[5]*g[1]*dg[k2]*dg[6 + k1] + w*d[7]*g[1]*dg[k2]*dg[6 + k1] + w*d[11]*g[1]*dg[k2]*dg[6 + k1] + w*d[15]*g[1]*dg[k2]*dg[6 + k1] + w*d[19]*g[1]*dg[k2]*dg[6 + k1] + w*d[21]*g[1]*dg[k2]*dg[6 + k1] + 2*w*d[8]*g[2]*dg[k2]*dg[6 + k1] + 2*w*d[20]*g[2]*dg[k2]*dg[6 + k1] + 2*w*d[24]*g[2]*dg[k2]*dg[6 + k1] + w*c[5]*dg[3 + k2]*dg[6 + k1] + w*c[7]*dg[3 + k2]*dg[6 + k1] + w*d[5]*g[0]*dg[3 + k2]*dg[6 + k1] + w*d[7]*g[0]*dg[3 + k2]*dg[6 + k1] + w*d[11]*g[0]*dg[3 + k2]*dg[6 + k1] + w*d[15]*g[0]*dg[3 + k2]*dg[6 + k1] + w*d[19]*g[0]*dg[3 + k2]*dg[6 + k1] + w*d[21]*g[0]*dg[3 + k2]*dg[6 + k1] + 2*w*d[14]*g[1]*dg[3 + k2]*dg[6 + k1] + 2*w*d[16]*g[1]*dg[3 + k2]*dg[6 + k1] + 2*w*d[22]*g[1]*dg[3 + k2]*dg[6 + k1] + 2*w*d[17]*g[2]*dg[3 + k2]*dg[6 + k1] + 2*w*d[23]*g[2]*dg[3 + k2]*dg[6 + k1] + 2*w*d[25]*g[2]*dg[3 + k2]*dg[6 + k1] + (dw[k1]*(b[2] + c[2]*g[0] + c[6]*g[0] + d[2]*std::pow(g[0],2) + d[6]*std::pow(g[0],2) + d[18]*std::pow(g[0],2) + c[5]*g[1] + c[7]*g[1] + d[5]*g[0]*g[1] + d[7]*g[0]*g[1] + d[11]*g[0]*g[1] + d[15]*g[0]*g[1] + d[19]*g[0]*g[1] + d[21]*g[0]*g[1] + d[14]*std::pow(g[1],2) + d[16]*std::pow(g[1],2) + d[22]*std::pow(g[1],2) + 2*(c[8] + (d[8] + d[20] + d[24])*g[0] + (d[17] + d[23] + d[25])*g[1])*g[2] + 3*d[26]*std::pow(g[2],2)) + w*(db[6 + k1] + g[0]*dc[6 + k1] + g[0]*dc[18 + k1] + g[1]*(dc[15 + k1] + dc[21 + k1]) + 2*g[2]*dc[24 + k1] + std::pow(g[0],2)*dd[6 + k1] + g[0]*g[1]*dd[15 + k1] + std::pow(g[0],2)*dd[18 + k1] + g[0]*g[1]*dd[21 + k1] + 2*g[0]*g[2]*dd[24 + k1] + g[0]*g[1]*dd[33 + k1] + std::pow(g[1],2)*dd[42 + k1] + g[0]*g[1]*dd[45 + k1] + std::pow(g[1],2)*dd[48 + k1] + 2*g[1]*g[2]*dd[51 + k1] + std::pow(g[0],2)*dd[54 + k1] + g[0]*g[1]*dd[57 + k1] + 2*g[0]*g[2]*dd[60 + k1] + g[0]*g[1]*dd[63 + k1] + std::pow(g[1],2)*dd[66 + k1] + 2*g[1]*g[2]*dd[69 + k1] + 2*g[0]*g[2]*dd[72 + k1] + 2*g[1]*g[2]*dd[75 + k1] + 3*std::pow(g[2],2)*dd[78 + k1] + (c[2] + c[6])*dg[k1] + 2*d[2]*g[0]*dg[k1] + 2*d[6]*g[0]*dg[k1] + 2*d[18]*g[0]*dg[k1] + d[5]*g[1]*dg[k1] + d[7]*g[1]*dg[k1] + d[11]*g[1]*dg[k1] + d[15]*g[1]*dg[k1] + d[19]*g[1]*dg[k1] + d[21]*g[1]*dg[k1] + 2*d[8]*g[2]*dg[k1] + 2*d[20]*g[2]*dg[k1] + 2*d[24]*g[2]*dg[k1] + (c[5] + c[7])*dg[3 + k1] + d[5]*g[0]*dg[3 + k1] + d[7]*g[0]*dg[3 + k1] + d[11]*g[0]*dg[3 + k1] + d[15]*g[0]*dg[3 + k1] + d[19]*g[0]*dg[3 + k1] + d[21]*g[0]*dg[3 + k1] + 2*d[14]*g[1]*dg[3 + k1] + 2*d[16]*g[1]*dg[3 + k1] + 2*d[22]*g[1]*dg[3 + k1] + 2*d[17]*g[2]*dg[3 + k1] + 2*d[23]*g[2]*dg[3 + k1] + 2*d[25]*g[2]*dg[3 + k1] + 2*c[8]*dg[6 + k1] + 2*((d[8] + d[20] + d[24])*g[0] + (d[17] + d[23] + d[25])*g[1] + 3*d[26]*g[2])*dg[6 + k1]))*dg[6 + k2];
      }
    }
    break;
  }
  return SymTensor(hess);
}

} // end namespace Spheral
