//---------------------------------Spheral++----------------------------------//
// CRKSPHUtilities
//
// Useful methods for using the CRKSPH formalism.
//
// Created by JMO, Fri Aug  8 16:16:33 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_CRKSPHFluidGradient__
#define __Spheral_NodeSpace_CRKSPHFluidGradient__

#include "CRKSPHCorrectionParams.hh"

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
             const CRKOrder correctionOrder,
             const typename Dimension::Vector& rij,
             const typename Dimension::Vector& etai,
             const typename Dimension::Scalar Hdeti,
             const typename Dimension::Vector& etaj,
             const typename Dimension::Scalar Hdetj,
             const typename Dimension::Scalar Ai,
             const typename Dimension::Vector& Bi,
             const typename Dimension::Tensor& Ci,
             const typename Dimension::Scalar correctionMin = std::numeric_limits<typename Dimension::Scalar>::lowest(),
             const typename Dimension::Scalar correctionMax = std::numeric_limits<typename Dimension::Scalar>::max());

// Compute the corrected kernel value, uncorrected and corrected gradients.
// Returned as the last three arguments.
template<typename Dimension>
void
CRKSPHKernelAndGradient(typename Dimension::Scalar& WCRKSPH,
                        typename Dimension::Scalar& gradWSPH,
                        typename Dimension::Vector& gradWCRKSPH,
                        const TableKernel<Dimension>& W,
                        const CRKOrder correctionOrder,
                        const typename Dimension::Vector& rij,
                        const typename Dimension::Vector& etai,
                        const typename Dimension::SymTensor& Hi,
                        const typename Dimension::Scalar Hdeti,
                        const typename Dimension::Vector& etaj,
                        const typename Dimension::SymTensor& Hj,
                        const typename Dimension::Scalar Hdetj,
                        const typename Dimension::Scalar Ai,
                        const typename Dimension::Vector& Bi,
                        const typename Dimension::Tensor& Ci,
                        const typename Dimension::Vector& gradAi,
                        const typename Dimension::Tensor& gradBi,
                        const typename Dimension::ThirdRankTensor& gradCi,
                        const typename Dimension::Scalar correctionMin = std::numeric_limits<typename Dimension::Scalar>::lowest(),
                        const typename Dimension::Scalar correctionMax = std::numeric_limits<typename Dimension::Scalar>::max());

}

#include "CRKSPHUtilitiesInline.hh"

#endif
