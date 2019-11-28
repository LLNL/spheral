//---------------------------------Spheral++----------------------------------//
// SuperiorRKUtilities
//
// Computes and evaluates RK corrections
//----------------------------------------------------------------------------//

#include "SuperiorRKUtilities.hh"

#include "Utilities/safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Evaluate the base kernel value, gradient, or hessian
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::Scalar
SuperiorRKUtilities<Dimension, correctionOrder>::
evaluateBaseKernel(const TableKernel<Dimension>& kernel,
                   const Vector& x,
                   const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto Hdet = H.Determinant();
  return kernel.kernelValue(etaMag, Hdet);
}

template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::Vector
SuperiorRKUtilities<Dimension, correctionOrder>::
evaluateBaseGradient(const TableKernel<Dimension>& kernel,
                     const Vector& x,
                     const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto etaMagInv = safeInv(etaMag);
  const auto Hdet = H.Determinant();
  const auto dk = kernel.gradValue(etaMag, Hdet);
  const auto HetaUnit = H * eta * etaMagInv;
  return HetaUnit * dk;
}

template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::SymTensor
SuperiorRKUtilities<Dimension, correctionOrder>::
evaluateBaseHessian(const TableKernel<Dimension>& kernel,
                    const Vector& x,
                    const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto etaMagInv = safeInv(etaMag);
  const auto Hdet = H.Determinant();
  const auto dk = kernel.gradValue(etaMag, Hdet);
  const auto ddk = kernel.grad2Value(etaMag, Hdet);
  const auto HetaUnit = H * eta * etaMagInv;
  const auto HetaUnit2 = HetaUnit.selfdyad();
  const auto H2 = H.square();
  return (H2 - HetaUnit2) * etaMagInv * dk + HetaUnit2 * ddk;
}

//------------------------------------------------------------------------------
// Evaluate the RK kernel, gradient, or Hessian
// W^{R}=C^{\top}PW
// W_{\gamma}^{R}=\left(C^{\top}P_{\gamma}+C_{\gamma}^{\top}P\right)W+C^{\top}PW_{\gamma}
// W_{\gamma\zeta}^{R}=\left(C^{\top}P_{\gamma\zeta}+C_{\gamma}^{\top}P_{\zeta}+C_{\zeta}^{\top}P_{\gamma}+C_{\gamma\zeta}^{\top}P\right)W+\left(C^{\top}P_{\gamma}+C_{\gamma}^{\top}P\right)W_{\zeta}+\left(C^{\top}P_{\zeta}+C_{\zeta}^{\top}P\right)W_{\gamma}+C^{\top}PW_{\gamma\zeta}
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::Scalar
SuperiorRKUtilities<Dimension, correctionOrder>::
evaluateKernel(const TableKernel<Dimension>& kernel,
               const Vector& x,
               const SymTensor& H,
               const std::vector<double>& corrections) {
  // Get kernel and polynomials
  const auto w = evaluateBaseKernel(kernel, x, H);
  const auto P = getPolynomials(x);
  
  // Get result
  const auto CP = innerProductRK(p, corrections, 0, 0);
  return CP * w;
}

template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::Vector
SuperiorRKUtilities<Dimension, correctionOrder>::
evaluateGradient(const TableKernel<Dimension>& kernel,
                 const Vector& x,
                 const SymTensor& H,
                 const std::vector<double>& corrections) {
  const auto dim = Dimension::nDim;
  
  // Get kernel and polynomials
  const auto w = evaluateBaseKernel(kernel, x, H);
  const auto dw = evaluteBaseGradient(kernel, x, H);
  const auto P = getPolynomials(x);
  const auto dP = getGradPolynomials(x);
  
  // Get result
  const auto CP = innerProductRK(p, corrections, 0, 0);
  Vector result = Vector::zero;
  for (auto d = 0; d < dimension; ++ d) {
    const auto CdP = innerProductRK(corrections, dP, 0, offsetGradP(d));
    const auto dCP = innerProductRK(corrections, P, offsetGradC(d), 0);
    result(d) = (Cdp + dCP) * w + CP * dw(d);
  }
  return result;
}

template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::SymTensor
SuperiorRKUtilities<Dimension, correctionOrder>::
evaluateHessian(const TableKernel<Dimension>& kernel,
                const Vector& x,
                const SymTensor& H,
                const std::vector<double>& corrections) {
  // Get kernel and polynomials
  const auto w = evaluateBaseKernel(kernel, x, H);
  const auto dw = evaluteBaseGradient(kernel, x, H);
  const auto ddw = evaluteBaseHessian(kernel, x, H);
  const auto P = getPolynomials(x);
  const auto dP = getGradPolynomials(x);
  const auto ddP = getHessPolynomials(x);

  // Get result
  const auto CP = innerProductRK(coefficients, P, 0, 0);
  SymTensor result = SymTensor::zero;
  for (auto d1 = 0; d1 < dimension; ++d1) {
    const auto Cd1P = innerProductRK(coefficients, dP, 0, offsetGradP(d1));
    const auto d1CP = innerProductRK(coefficients, P, offsetGradC(d1), 0);
    for (auto d2 = d1; d2 < dimension; ++d2) {
      const auto Cd2P = innerProductRK(coefficients, ddP, 0, offsetGradP(d2));
      const auto d2CP = innerProductRK(coefficients, P, offsetGradC(d2), 0);
      const auto CddP = innerProductRK(coefficients, ddP, 0, offsetHessP(d1, d2));
      const auto d1Cd2P = innerProductRK(coefficients, dP, offsetGradC(d1), offsetGradP(d2));
      const auto d2Cd1P = innerProductRK(coefficients, dP, offsetGradC(d2), offsetGradP(d1));
      const auto ddCP = innerProductRK(coefficients, P, offsetHessC(d1, d2), 0);
      
      result(d1, d2) = (CddP + d1Cd2P + d2Cd1P + ddCP) * w + (Cd1P + d1CP) * dw(d2) + (Cd2P + d2CP) * dw(d1) + CP * dw(d1, d2);
    }
  }
  return result;
}



} // end namespace Spheral
