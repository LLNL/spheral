//---------------------------------Spheral++----------------------------------//
// CRKSPHUtilities
//
// Useful methods for using the CRKSPH formalism.
//
// Created by JMO, Fri Aug  8 16:16:33 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_CRKSPHFluidGradient__
#define __Spheral_NodeSpace_CRKSPHFluidGradient__

#include "RK/RKCorrectionParams.hh"

#include <limits>

// Forward declarations.
namespace Spheral {
  template<typename Dimension> class TableKernel;
}

namespace Spheral {

// Compute the corrected kernel value.
template<typename Dimension>
typename Dimension::Scalar
CRKSPHKernel(const TableKernel<Dimension>& W,
             const RKOrder correctionOrder,
             const typename Dimension::Vector& rij,
             const typename Dimension::Vector& etaj,
             const typename Dimension::Scalar Hdetj,
             const typename Dimension::Scalar Ai,
             const typename Dimension::Vector& Bi,
             const typename Dimension::Tensor& Ci);

// Compute the corrected kernel value, uncorrected and corrected gradients.
// Returned as the last three arguments.
template<typename Dimension>
void
CRKSPHKernelAndGradient(typename Dimension::Scalar& WCRKSPH,
                        typename Dimension::Scalar& gradWSPH,
                        typename Dimension::Vector& gradWCRKSPH,
                        const TableKernel<Dimension>& W,
                        const RKOrder correctionOrder,
                        const typename Dimension::Vector& rij,
                        const typename Dimension::Vector& etaj,
                        const typename Dimension::SymTensor& Hj,
                        const typename Dimension::Scalar Hdetj,
                        const typename Dimension::Scalar Ai,
                        const typename Dimension::Vector& Bi,
                        const typename Dimension::Tensor& Ci,
                        const typename Dimension::Vector& gradAi,
                        const typename Dimension::Tensor& gradBi,
                        const typename Dimension::ThirdRankTensor& gradCi);

}

#include "CRKSPHUtilitiesInline.hh"

#endif
