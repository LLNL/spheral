//---------------------------------Spheral++----------------------------------//
// CRKSPHUtilities
//
// Useful methods for using the CRKSPH formalism.
//
// Created by JMO, Fri Aug  8 16:16:33 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_CRKSPHFluidGradient__
#define __Spheral_NodeSpace_CRKSPHFluidGradient__

// Forward declarations.
namespace Spheral {
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
}

namespace Spheral {
namespace CRKSPHSpace {

// Compute the corrected kernel value.
template<typename Dimension>
typename Dimension::Scalar
CRKSPHKernel(const KernelSpace::TableKernel<Dimension>& W,
           const typename Dimension::Vector& rij,
           const typename Dimension::Vector& etai,
           const typename Dimension::Scalar& Hdeti,
           const typename Dimension::Vector& etaj,
           const typename Dimension::Scalar& Hdetj,
           const typename Dimension::Scalar& Ai,
           const typename Dimension::Vector& Bi);

// Compute the corrected kernel value, uncorrected and corrected gradients.
// Returned as the last three arguments.
template<typename Dimension>
void
CRKSPHKernelAndGradient(const KernelSpace::TableKernel<Dimension>& W,
                      const typename Dimension::Vector& rij,
                      const typename Dimension::Vector& etai,
                      const typename Dimension::SymTensor& Hi,
                      const typename Dimension::Scalar& Hdeti,
                      const typename Dimension::Vector& etaj,
                      const typename Dimension::SymTensor& Hj,
                      const typename Dimension::Scalar& Hdetj,
                      const typename Dimension::Scalar& Ai,
                      const typename Dimension::Vector& Bi,
                      const typename Dimension::Vector& gradAi,
                      const typename Dimension::Tensor& gradBi,
                      typename Dimension::Scalar& WCRKSPH,
                      typename Dimension::Scalar& gradWSPH,
                      typename Dimension::Vector& gradWCRKSPH);

}
}

#ifndef __GCCXML__
#include "CRKSPHUtilitiesInline.hh"
#endif

#endif
