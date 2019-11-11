//---------------------------------Spheral++----------------------------------//
// CRKSPHUtilities
//
// Useful methods for using the CRKSPH formalism.
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_CRKSPHFluidGradient__
#define __Spheral_NodeSpace_CRKSPHFluidGradient__

#include "CRKSPHCorrectionParams.hh"

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class TableKernel;
}

namespace Spheral {

// Compute the corrected kernel value.
template<typename Dimension>
typename Dimension::Scalar
evaluateRKKernel(const TableKernel<Dimension>& W,
                 const CRKOrder correctionOrder,
                 const typename Dimension::Vector& rij,
                 const typename Dimension::Vector& etaj,
                 const typename Dimension::Scalar Hdetj,
                 const typename Dimension::Scalar Ai,
                 const typename Dimension::Vector& Bi,
                 const typename Dimension::Tensor& Ci,
                 const typename Dimension::ThirdRankTensor& Di);

// Compute the corrected kernel gradient.
template<typename Dimension>
typename Dimension::Vector
evaluateRKGradient(const TableKernel<Dimension>& W,
                   const CRKOrder correctionOrder,
                   const typename Dimension::Vector& rij,
                   const typename Dimension::Vector& etaj,
                   const typename Dimension::SymTensor& Hj,
                   const typename Dimension::Scalar Hdetj,
                   const typename Dimension::Scalar Ai,
                   const typename Dimension::Vector& Bi,
                   const typename Dimension::Tensor& Ci,
                   const typename Dimension::ThirdRankTensor& Di,
                   const typename Dimension::Vector& gradAi,
                   const typename Dimension::Tensor& gradBi,
                   const typename Dimension::ThirdRankTensor& gradCi,
                   const typename Dimension::FourthRankTensor& gradDi);

// Compute the corrected kernel hessian.
template<typename Dimension>
typename Dimension::Tensor
evaluateRKHessian(const TableKernel<Dimension>& W,
                  const CRKOrder correctionOrder,
                  const typename Dimension::Vector& etaj,
                  const typename Dimension::SymTensor& Hj,
                  const typename Dimension::Scalar Hdetj,
                  const typename Dimension::Scalar Ai,
                  const typename Dimension::Vector& Bi,
                  const typename Dimension::Tensor& Ci,
                  const typename Dimension::ThirdRankTensor& Di,
                  const typename Dimension::Vector& gradAi,
                  const typename Dimension::Tensor& gradBi,
                  const typename Dimension::ThirdRankTensor& gradCi,
                  const typename Dimension::FourthRankTensor& gradDi,
                  const typename Dimension::Tensor& hessAi,
                  const typename Dimension::ThirdRankTensor& hessBi,
                  const typename Dimension::FourthRankTensor& hessCi,
                  const typename Dimension::FifthRankTensor& hessDi);
}

#include "CRKSPHUtilitiesInline.hh"

#endif
