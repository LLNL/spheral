//---------------------------------Spheral++----------------------------------//
// CRKSPHUtilities
//
// Useful methods for using the CRKSPH formalism.
//----------------------------------------------------------------------------//
#include "Kernel/TableKernel.hh"
#include "CRKSPHCorrectionParams.hh"
#include "Geometry/Dimension.hh"
#include "safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the corrected kernel value
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
evaluateRKKernel(const TableKernel<Dimension>& kernel,
                 const CRKOrder correctionOrder,
                 const typename Dimension::Vector& eta,
                 const typename Dimension::SymTensor& H,
                 const typename Dimension::Scalar a,
                 const typename Dimension::Vector& b,
                 const typename Dimension::Tensor& c,
                 const typename Dimension::ThirdRankTensor& d) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

  // Evaluate the kernel
  const auto w = kernel(eta, H);

  // Initialize values
  Scalar val = Scalar::zero;
  constexpr auto dim = Dimension::nDim;

  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::Cubic:
    for (auto q1 = 0; q1 < dim; ++q1) {
      for (auto q2 = 0; q2 < dim; ++q2) {
        for (auto q3 = 0; q3 < dim; ++q3) {
          val += d[q1][q2][q3]*eta[q1]*eta[q2]*eta[q3];
        }
      }
    }
  case CRKOrder::Quadratic:
    for (auto q1 = 0; q1 < dim; ++q1) {
      for (auto q2 = 0; q2 < dim; ++q2) {
        val += c[q1][q2]*eta[q1]*eta[q2];
      }
    }
  case CRKOrder::ZerothOrder:
    for (auto q1 = 0; q1 < dim; ++q1) {
      val += b[q1]*eta[q1];
    }
  case CRKOrder::Linear:
    val += a;
    break;
  }
  
  return val * w;
}

//------------------------------------------------------------------------------
// Compute the corrected kernel gradient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
evaluateRKGradient(const TableKernel<Dimension>& kernel,
                   const CRKOrder correctionOrder,
                   const typename Dimension::Vector& eta,
                   const typename Dimension::SymTensor& H,
                   const typename Dimension::Scalar a,
                   const typename Dimension::Vector& b,
                   const typename Dimension::Tensor& c,
                   const typename Dimension::ThirdRankTensor& d,
                   const typename Dimension::Vector& da,
                   const typename Dimension::Tensor& db,
                   const typename Dimension::ThirdRankTensor& dc,
                   const typename Dimension::FourthRankTensor& dd) {
  typedef typename Dimension::Vector Vector;

  // Evaluate the kernel
  const auto deta = H;
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto c = kernel(eta, H);
  const auto dc = kernel.grad(eta, H);
  const auto ddc = kernel.grad2(eta, H);
  const auto w = c;
  const auto dw = Heta * dc;
  
  // Initialize values
  Vector grad = Vector::zero;
  constexpr auto dim = Dimension::nDim;
  
  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::Cubic:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        for (auto q2 = 0; q2 < dim; ++q2) {
          for (auto q3 = 0; q3 < dim; ++q3) {
            grad[k1] += dw[k1]*(d[q1][q2][q3]*eta[q1]*eta[q2]*eta[q3]) + w*(eta[q1]*eta[q2]*eta[q3]*dd[q1][q2][q3][k1] + d[q1][q2][q3]*eta[q2]*eta[q3]*deta[q1][k1] + d[q1][q2][q3]*eta[q1]*eta[q3]*deta[q2][k1] + d[q1][q2][q3]*eta[q1]*eta[q2]*deta[q3][k1]);
          }
        }
      }
    }
  case CRKOrder::Quadratic:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        for (auto q2 = 0; q2 < dim; ++q2) {
          grad[k1] += dw[k1]*c[q1][q2]*eta[q1]*eta[q2] + w*(eta[q1]*eta[q2]*dc[q1][q2][k1] + c[q1][q2]*eta[q2]*deta[q1][k1] + c[q1][q2]*eta[q1]*deta[q2][k1]);
        }
      }
    }
  case CRKOrder::Linear:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto q1 = 0; q1 < dim; ++q1) {
        grad[k1] += dw[k1]*b[q1]*eta[q1] + w*(eta[q1]*db[q1][k1] + b[q1]*deta[q1][k1]);
      }
    }
  case CRKOrder::ZerothOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      grad[k1] += dw[k1]*a + w*da[k1];
    }
  }

  return grad;
}

//------------------------------------------------------------------------------
// Compute the corrected kernel gradient
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
evaluateRKHessian(const TableKernel<Dimension>& kernel,
                  const CRKOrder correctionOrder,
                  const typename Dimension::Vector& eta,
                  const typename Dimension::SymTensor& H,
                  const typename Dimension::Scalar a,
                  const typename Dimension::Vector& b,
                  const typename Dimension::Tensor& c,
                  const typename Dimension::ThirdRankTensor& d,
                  const typename Dimension::Vector& da,
                  const typename Dimension::Tensor& db,
                  const typename Dimension::ThirdRankTensor& dc,
                  const typename Dimension::FourthRankTensor& dd,
                  const typename Dimension::Tensor& dda,
                  const typename Dimension::ThirdRankTensor& ddb,
                  const typename Dimension::FourthRankTensor& ddc,
                  const typename Dimension::FifthRankTensor& ddd) {
  typedef typename Dimension::Tensor Tensor;
  typedef typename Dimension::ThirdRankTensor ThirdRankTensor;

  // Evaluate the kernel
  const auto deta = H;
  const auto ddeta = ThirdRankTensor::zero;
  const auto etaMagInv = safeInv(eta.magnitude());
  const auto Heta = H * eta * etaMagInv;
  const auto Heta2 = Heta.selfdyad();
  const auto H2 = H.square();
  const auto c = kernel(eta, H);
  const auto dc = kernel.grad(eta, H);
  const auto ddc = kernel.grad2(eta, H);
  const auto w = c;
  const auto dw = Heta * dc;
  const auto ddw = (H2 - Heta2) * etaMagInv * dc + Heta2 * ddc;
  
  // Initialize values
  Tensor hess = Tensor::zero;
  constexpr auto dim = Dimension::nDim;

  // Return kernel value
  switch (correctionOrder) {
  case CRKOrder::Cubic:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        for (auto q1 = 0; q1 < dim; ++q1) {
          for (auto q2 = 0; q2 < dim; ++q2) {
            for (auto q3 = 0; q3 < dim; ++q3) {
              hess[k1][k2] += ddw[k1][k2]*d[q1][q2][q3]*eta[q1]*eta[q2]*eta[q3] + dw[k2]*(eta[q1]*eta[q2]*eta[q3]*dd[q1][q2][q3][k1] + d[q1][q2][q3]*eta[q2]*eta[q3]*deta[q1][k1] + d[q1][q2][q3]*eta[q1]*eta[q3]*deta[q2][k1] + d[q1][q2][q3]*eta[q1]*eta[q2]*deta[q3][k1]) + dw[k1]*(eta[q1]*eta[q2]*eta[q3]*dd[q1][q2][q3][k2] + d[q1][q2][q3]*eta[q2]*eta[q3]*deta[q1][k2] + d[q1][q2][q3]*eta[q1]*eta[q3]*deta[q2][k2] + d[q1][q2][q3]*eta[q1]*eta[q2]*deta[q3][k2]) + w*(eta[q1]*eta[q2]*eta[q3]*ddd[q1][q2][q3][k1][k2] + d[q1][q2][q3]*eta[q2]*eta[q3]*ddeta[q1][k1][k2] + d[q1][q2][q3]*eta[q1]*eta[q3]*ddeta[q2][k1][k2] + d[q1][q2][q3]*eta[q1]*eta[q2]*ddeta[q3][k1][k2] + eta[q2]*eta[q3]*dd[q1][q2][q3][k2]*deta[q1][k1] + eta[q2]*eta[q3]*dd[q1][q2][q3][k1]*deta[q1][k2] + eta[q1]*eta[q3]*dd[q1][q2][q3][k2]*deta[q2][k1] + d[q1][q2][q3]*eta[q3]*deta[q1][k2]*deta[q2][k1] + eta[q1]*eta[q3]*dd[q1][q2][q3][k1]*deta[q2][k2] + d[q1][q2][q3]*eta[q3]*deta[q1][k1]*deta[q2][k2] + eta[q1]*eta[q2]*dd[q1][q2][q3][k2]*deta[q3][k1] + d[q1][q2][q3]*eta[q2]*deta[q1][k2]*deta[q3][k1] + d[q1][q2][q3]*eta[q1]*deta[q2][k2]*deta[q3][k1] + eta[q1]*eta[q2]*dd[q1][q2][q3][k1]*deta[q3][k2] + d[q1][q2][q3]*eta[q2]*deta[q1][k1]*deta[q3][k2] + d[q1][q2][q3]*eta[q1]*deta[q2][k1]*deta[q3][k2]);
            }
          }
        }
      }
    }
  case CRKOrder::Quadratic:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        for (auto q1 = 0; q1 < dim; ++q1) {
          for (auto q2 = 0; q2 < dim; ++q2) {
            hess[k1][k2] += ddw[k1][k2]*c[q1][q2]*eta[q1]*eta[q2] + dw[k2]*(eta[q1]*eta[q2]*dc[q1][q2][k1] + c[q1][q2]*eta[q2]*deta[q1][k1] + c[q1][q2]*eta[q1]*deta[q2][k1]) + dw[k1]*(da[k2] + eta[q1]*eta[q2]*dc[q1][q2][k2] + c[q1][q2]*eta[q2]*deta[q1][k2] + c[q1][q2]*eta[q1]*deta[q2][k2]) + w*(eta[q1]*eta[q2]*ddc[q1][q2][k1][k2] + c[q1][q2]*eta[q2]*ddeta[q1][k1][k2] + c[q1][q2]*eta[q1]*ddeta[q2][k1][k2] + eta[q2]*dc[q1][q2][k2]*deta[q1][k1] + eta[q2]*dc[q1][q2][k1]*deta[q1][k2] + eta[q1]*dc[q1][q2][k2]*deta[q2][k1] + c[q1][q2]*deta[q1][k2]*deta[q2][k1] + eta[q1]*dc[q1][q2][k1]*deta[q2][k2] + c[q1][q2]*deta[q1][k1]*deta[q2][k2]);
          }
        }
      }
    }
  case CRKOrder::Linear:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        for (auto q1 = 0; q1 < dim; ++q1) {
          hess[k1][k2] += ddw[k1][k2]*b[q1]*eta[q1] + dw[k2]*(eta[q1]*db[q1][k1] + b[q1]*deta[q1][k1]) + dw[k1]*(eta[q1]*db[q1][k2] + b[q1]*deta[q1][k2]) + w*(eta[q1]*ddb[q1][k1][k2] + b[q1]*ddeta[q1][k1][k2] + db[q1][k2]*deta[q1][k1] + db[q1][k1]*deta[q1][k2]);
        }
      }
    }
  case CRKOrder::ZerothOrder:
    for (auto k1 = 0; k1 < dim; ++k1) {
      for (auto k2 = 0; k2 < dim; ++k2) {
        hess[k1][k2] += w*dda[k1][k2] + a*ddw[k1][k2] + da[k2]*dw[k1] + da[k1]*dw[k2]
      }
    }
  }

  return hess;
}

} // end namespace Spheral
