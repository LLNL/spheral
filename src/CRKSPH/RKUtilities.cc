//---------------------------------Spheral++----------------------------------//
// RKUtilities
//
// Computes and evaluates RK corrections
//----------------------------------------------------------------------------//

#include "RKUtilities.hh"
#include "Eigen/Dense"
#include "Neighbor/ConnectivityMap.hh"
#include "Utilities/safeInv.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Evaluate the base kernel value, gradient, or hessian
//------------------------------------------------------------------------------
template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::Scalar
RKUtilities<Dimension, correctionOrder>::
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
RKUtilities<Dimension, correctionOrder>::
evaluateBaseGradient(const TableKernel<Dimension>& kernel,
                     const Vector& x,
                     const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto etaUnit = eta.unitVector();
  const auto Hdet = H.Determinant();
  const auto dk = kernel.gradValue(etaMag, Hdet);
  const auto HetaUnit = H * etaUnit;
  return HetaUnit * dk;
}

template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::SymTensor
RKUtilities<Dimension, correctionOrder>::
evaluateBaseHessian(const TableKernel<Dimension>& kernel,
                    const Vector& x,
                    const SymTensor& H) {
  const auto eta = H * x;
  const auto etaMag = eta.magnitude();
  const auto etaUnit = eta.unitVector();
  const auto etaMagInv = safeInv(etaMag);
  const auto Hdet = H.Determinant();
  const auto dk = kernel.gradValue(etaMag, Hdet);
  const auto ddk = kernel.grad2Value(etaMag, Hdet);
  const auto HetaUnit = H * etaUnit;
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
RKUtilities<Dimension, correctionOrder>::
evaluateKernel(const TableKernel<Dimension>& kernel,
               const Vector& x,
               const SymTensor& H,
               const std::vector<double>& corrections) {
  CHECK(corrections.size() == correctionsSize(false)
        || corrections.size() == correctionsSize(true));
  
  // Get kernel and polynomials
  const auto w = evaluateBaseKernel(kernel, x, H);
  const auto P = getPolynomials(x);
  
  // Perform inner products and return result
  const auto CP = innerProductRK(P, corrections, 0, 0);
  return CP * w;
}

template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::Vector
RKUtilities<Dimension, correctionOrder>::
evaluateGradient(const TableKernel<Dimension>& kernel,
                 const Vector& x,
                 const SymTensor& H,
                 const std::vector<double>& corrections) {
  CHECK(corrections.size() == correctionsSize(false)
        || corrections.size() == correctionsSize(true));
  
  const auto dim = Dimension::nDim;
  
  // Get kernel and polynomials
  const auto w = evaluateBaseKernel(kernel, x, H);
  const auto dw = evaluateBaseGradient(kernel, x, H);
  const auto P = getPolynomials(x);
  const auto dP = getGradPolynomials(x);
  
  // Perform inner products and return result
  const auto CP = innerProductRK(corrections, P, 0, 0);
  Vector result = Vector::zero;
  for (auto d = 0; d < dim; ++ d) {
    const auto CdP = innerProductRK(corrections, dP, 0, offsetGradP(d));
    const auto dCP = innerProductRK(corrections, P, offsetGradC(d), 0);
    result(d) = (CdP + dCP) * w + CP * dw(d);
  }
  return result;
}

template<typename Dimension, CRKOrder correctionOrder>
typename Dimension::SymTensor
RKUtilities<Dimension, correctionOrder>::
evaluateHessian(const TableKernel<Dimension>& kernel,
                const Vector& x,
                const SymTensor& H,
                const std::vector<double>& corrections) {
  CHECK(corrections.size() == correctionsSize(false)
        || corrections.size() == correctionsSize(true));
  
  const auto dim = Dimension::nDim;
  
  // Get kernel and polynomials
  const auto w = evaluateBaseKernel(kernel, x, H);
  const auto dw = evaluateBaseGradient(kernel, x, H);
  const auto ddw = evaluateBaseHessian(kernel, x, H);
  const auto P = getPolynomials(x);
  const auto dP = getGradPolynomials(x);
  const auto ddP = getHessPolynomials(x);

  // Perform inner products and return result
  // Could precompute the inner products in a separate loop for efficiency
  const auto CP = innerProductRK(corrections, P, 0, 0);
  SymTensor result = SymTensor::zero;
  for (auto d1 = 0; d1 < dim; ++d1) {
    const auto Cd1P = innerProductRK(corrections, dP, 0, offsetGradP(d1));
    const auto d1CP = innerProductRK(corrections, P, offsetGradC(d1), 0);
    for (auto d2 = d1; d2 < dim; ++d2) {
      const auto Cd2P = innerProductRK(corrections, ddP, 0, offsetGradP(d2));
      const auto d2CP = innerProductRK(corrections, P, offsetGradC(d2), 0);
      const auto CddP = innerProductRK(corrections, ddP, 0, offsetHessP(d1, d2));
      const auto d1Cd2P = innerProductRK(corrections, dP, offsetGradC(d1), offsetGradP(d2));
      const auto d2Cd1P = innerProductRK(corrections, dP, offsetGradC(d2), offsetGradP(d1));
      const auto ddCP = innerProductRK(corrections, P, offsetHessC(d1, d2), 0);
      
      result(d1, d2) = (CddP + d1Cd2P + d2Cd1P + ddCP) * w + (Cd1P + d1CP) * dw(d2) + (Cd2P + d2CP) * dw(d1) + CP * ddw(d1, d2);
    }
  }
  return result;
}

// Compute the corrections
template<typename Dimension, CRKOrder correctionOrder>
void
RKUtilities<Dimension, correctionOrder>::
computeCorrections(const ConnectivityMap<Dimension>& connectivityMap,
                   const TableKernel<Dimension>& kernel,
                   const FieldList<Dimension, Scalar>& volume,
                   const FieldList<Dimension, Vector>& position,
                   const FieldList<Dimension, SymTensor>& H,
                   const bool needHessian,
                   FieldList<Dimension, std::vector<double>>& corrections) {
  // Typedefs: Eigen requires aligned allocator for stl containers before c++17
  typedef Eigen::Matrix<double, polynomialSize, 1> VectorType;
  typedef Eigen::Matrix<double, polynomialSize, polynomialSize> MatrixType;
  typedef std::vector<VectorType, Eigen::aligned_allocator<VectorType>> VectorOfVectorType;
  typedef std::vector<MatrixType, Eigen::aligned_allocator<MatrixType>> VectorOfMatrixType;
  
  // Size info
  const auto dim = Dimension::nDim;
  const auto numNodeLists = volume.size();
  const auto hessSize = needHessian ? symmetricMatrixSize(dim) : 0;
  
  // Check things
  REQUIRE(position.size() == numNodeLists);
  REQUIRE(H.size() == numNodeLists);
  REQUIRE(corrections.size() == numNodeLists);

  // Add values to M for a given j
  auto addToM = [&dim](const Scalar& v,
                       const Scalar& w,
                       const std::vector<double>& p,
                       MatrixType& M) {
    for (auto k = 0; k < polynomialSize; ++k) {
      for (auto l = k; l < polynomialSize; ++l) {
        M(k, l) += v * p[k] * p[l] * w;
      }
    }
    return;
  };
  
  // Add values to gradM for a given j
  auto addTodM = [&dim](const Scalar& v,
                        const Scalar& w,
                        const Vector& dw,
                        const std::vector<double>& p,
                        const std::vector<double>& dp,
                        VectorOfMatrixType& dM) {
    for (auto d = 0; d < dim; ++d) {
      const auto offd = offsetGradP(d);
      for (auto k = 0; k < polynomialSize; ++k) {
        for (auto l = k; l < polynomialSize; ++l) {
          dM[d](k,l) += v * ((dp[offd+k] * p[l] + p[k] * dp[offd+l]) * w + p[k] * p[l] * dw(d));
        }
      }
    }
    return;
  };

  // Add values to hessM for a given j
  auto addToddM = [&dim](const Scalar& v,
                         const Scalar& w,
                         const Vector& dw,
                         const SymTensor& ddw,
                         const std::vector<double>& p,
                         const std::vector<double>& dp,
                         const std::vector<double>& ddp,
                         VectorOfMatrixType& ddM) {
    for (auto d1 = 0; d1 < dim; ++d1) {
      const auto offd1 = offsetGradP(d1);
      for (auto d2 = d1; d2 < dim; ++d2) {
        const auto offd2 = offsetGradP(d2);
        const auto offd12 = offsetHessP(d1, d2);
        const auto d12 = flatSymmetricIndex(d1, d2);
        for (auto k = 0; k < polynomialSize; ++k) {
          for (auto l = k; l < polynomialSize; ++l) {
            ddM[d12](k,l) += v * ((ddp[offd12+k] * p[l] + dp[offd1+k] * dp[offd2+l] + dp[offd2+k] * dp[offd1+l] + p[k] * ddp[offd12+l]) * w + (dp[offd1+k] * p[l] + p[k] * dp[offd1+l]) * dw(d2) + (dp[offd2+k] * p[l] + p[k] * dp[offd2+l]) * dw(d1) + p[k] * p[l] * ddw(d1, d2));
          }
        }
      }
    }
    return;
  };

  // Compute corrections for each point independently
  MatrixType M;
  VectorOfMatrixType dM(dim);
  VectorOfMatrixType ddM(hessSize);
  VectorType rhs;
  VectorType C;
  VectorOfVectorType dC(dim);
  VectorOfVectorType ddC(hessSize);
  for (auto nodeListi = 0; nodeListi < numNodeLists; ++nodeListi) {
    const auto numNodes = connectivityMap.numNodes(nodeListi);
    for (auto nodei = 0; nodei < numNodes; ++nodei) {
      // Get data for point i
      const auto xi = position(nodeListi , nodei);
      
      // Initialize polynomial matrices for point i
      M.setZero();
      for (auto& mat : dM) {
        mat.setZero();
      }
      for (auto& mat : ddM) {
        mat.setZero();
      }
      
      // Get function for adding contribution to matrices
      auto addToMatrix = [&](const int nodeListj,
                             const int nodej) {
        // Get data for point j
        const auto xj = position(nodeListj, nodej);
        const auto xij = xi - xj;
        const auto Hj = H(nodeListj, nodej);
        const auto vj = volume(nodeListj, nodej);
        
        // Add to matrix
        const auto w = evaluateBaseKernel(kernel, xij, Hj);
        const auto p = getPolynomials(xij);
        CHECK(p.size() == polynomialSize);
        addToM(vj, w, p, M);

        // Add to gradient matrix
        const auto dw = evaluateBaseGradient(kernel, xij, Hj);
        const auto dp = getGradPolynomials(xij);
        CHECK(dp.size() == polynomialSize * dim);
        addTodM(vj, w, dw, p, dp, dM);
        
        // Add to Hessian matrix
        if (needHessian) {
          const auto ddw = evaluateBaseHessian(kernel, xij, Hj);
          const auto ddp = getHessPolynomials(xij);
          CHECK(ddp.size() == hessSize * polynomialSize);
          addToddM(vj, w, dw, ddw, p, dp, ddp, ddM);
        }
        
        return;
      };
                            
      // Add contribution from other points
      const auto& connectivity = connectivityMap.connectivityForNode(nodeListi, nodei);
      for (auto nodeListj = 0; nodeListj < numNodeLists; ++nodeListj) {
        for (auto nodej : connectivity[nodeListj]) {
          addToMatrix(nodeListj, nodej);
        } // nodej
      } // nodeListj
      
      // Add self contribution
      addToMatrix(nodeListi, nodei);
      
      // M symmetries
      for (auto k = 0; k < polynomialSize; ++k) {
        for (auto l = 0; l < k; ++l) {
          M(k, l) = M(l, k);
        }
      }

      // dM symmetries
      for (auto d = 0; d < dim; ++d) {
        for (auto k = 0; k < polynomialSize; ++k) {
          for (auto l = 0; l < k; ++l) {
            dM[d](k, l) = dM[d](l, k);
          }
        }
      }

      // ddM symmetries
      if (needHessian) {
        for (auto d1 = 0; d1 < dim; ++d1) {
          for (auto d2 = 0; d2 < dim; ++d2) {
            const auto d12 = flatSymmetricIndex(d1, d2);
            for (auto k = 0; k < polynomialSize; ++k) {
              for (auto l = 0; l < k; ++l) {
                ddM[d12](k,l) = ddM[d12](l,k);
              }
            }
          }
        }
      }
      
      // Get inverse of M matrix
      auto solver = M.colPivHouseholderQr();
      
      // Compute corrections
      rhs(0) = 1;
      for (auto k = 1; k < polynomialSize; ++k) {
        rhs(k) = 0;
      }
      C = solver.solve(rhs);

      // Compute gradient corrections
      for (auto d = 0; d < dim; ++d) {
        rhs = -(dM[d] * C);
        dC[d] = solver.solve(rhs);
      }
      
      // Compute hessian corrections
      if (needHessian) {
        for (auto d1 = 0; d1 < dim; ++d1) {
          for (auto d2 = 1; d2 < dim; ++d2) {
            const auto d12 = flatSymmetricIndex(d1, d2);
            rhs = -(ddM[d12] * C + dM[d1] * dC[d2] + dM[d2] * dC[d1]);
            ddC[d12] = solver.solve(rhs);
          }
        }
      }

      // Initialize corrections vector
      const auto corrSize = correctionsSize(needHessian);
      auto& corr = corrections(nodeListi, nodei);
      corr.resize(corrSize);
      
      // Put corrections into vector
      for (auto k = 0; k < polynomialSize; ++k) {
        corr[k] = C(k);
      }

      // Put gradient corrections into vector
      for (auto d = 0; d < dim; ++d) {
        const auto offd = offsetGradC(d);
        for (auto k = 0; k < polynomialSize; ++k) {
          corr[offd+k] = dC[d](k);
        }
      }

      // Put hessian corrections into vector
      if (needHessian) {
        for (auto d1 = 0; d1 < dim; ++d1) {
          for (auto d2 = d1; d2 < dim; ++d2) {
            const auto d12 = flatSymmetricIndex(d1, d2);
            const auto offd12 = offsetHessC(d1, d2);
            for (auto k = 0; k < polynomialSize; ++k) {
              corr[offd12+k] = ddC[d12](k);
            }
          }
        }
      }
    } // nodei
  } // nodeListi
} // computeCorrections

} // end namespace Spheral
