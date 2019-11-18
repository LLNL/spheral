//---------------------------------Spheral++----------------------------------//
// RKUtilities
//
// Evaluate the RK kernel
//----------------------------------------------------------------------------//
#ifndef __Spheral_RKUtilities__
#define __Spheral_RKUtilities__

#include "CRKSPHCorrectionParams.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class TableKernel;
}

namespace Spheral {

// Compute the corrected kernel value.
template<typename Dimension>
typename Dimension::Scalar
evaluateRKKernel(const TableKernel<Dimension>& kernel,
                 const CRKOrder correctionOrder,
                 const typename Dimension::Vector& eta,
                 const typename Dimension::SymTensor& H,
                 const typename Dimension::Scalar& a,
                 const typename Dimension::Vector& b,
                 const typename Dimension::Tensor& c,
                 const typename Dimension::ThirdRankTensor& d);

// Compute the corrected kernel gradient.
template<typename Dimension>
typename Dimension::Vector
evaluateRKGradient(const TableKernel<Dimension>& kernel,
                   const CRKOrder correctionOrder,
                   const typename Dimension::Vector& eta,
                   const typename Dimension::SymTensor& H,
                   const typename Dimension::Scalar& a,
                   const typename Dimension::Vector& b,
                   const typename Dimension::Tensor& c,
                   const typename Dimension::ThirdRankTensor& d,
                   const typename Dimension::Vector& da,
                   const typename Dimension::Tensor& db,
                   const typename Dimension::ThirdRankTensor& dc,
                   const typename Dimension::FourthRankTensor& dd);

// Compute the corrected kernel hessian.
template<typename Dimension>
typename Dimension::SymTensor
evaluateRKHessian(const TableKernel<Dimension>& kernel,
                  const CRKOrder correctionOrder,
                  const typename Dimension::Vector& eta,
                  const typename Dimension::SymTensor& H,
                  const typename Dimension::Scalar& a,
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
                  const typename Dimension::FifthRankTensor& ddd);

} // end namespace Spheral

#endif
