//---------------------------------Spheral++----------------------------------//
// CSPHUtilities
//
// Useful methods for using the CSPH formalism.
//
// Created by JMO, Fri Aug  8 16:16:33 PDT 2008
//----------------------------------------------------------------------------//
#ifndef __Spheral_NodeSpace_CSPHFluidGradient__
#define __Spheral_NodeSpace_CSPHFluidGradient__

// Forward declarations.
namespace Spheral {
  namespace KernelSpace {
    template<typename Dimension> class TableKernel;
  }
}

namespace Spheral {
namespace CSPHSpace {

// Compute the corrected kernel value.
template<typename Dimension>
typename Dimension::Scalar
CSPHKernel(const KernelSpace::TableKernel<Dimension>& W,
           const typename Dimension::Vector& rij,
           const typename Dimension::Vector& etaj,
           const typename Dimension::Scalar& Hdetj,
           const typename Dimension::Scalar& Ai,
           const typename Dimension::Vector& Bi);

// Compute the corrected kernel value, uncorrected and corrected gradients.
// Returned as the last three arguments.
template<typename Dimension>
void
CSPHKernelAndGradient(const KernelSpace::TableKernel<Dimension>& W,
                      const typename Dimension::Vector& rij,
                      const typename Dimension::Vector& etaj,
                      const typename Dimension::SymTensor& Hj,
                      const typename Dimension::Scalar& Hdetj,
                      const typename Dimension::Scalar& Ai,
                      const typename Dimension::Vector& Bi,
                      const typename Dimension::Vector& gradAi,
                      const typename Dimension::Tensor& gradBi,
                      typename Dimension::Scalar& WCSPH,
                      typename Dimension::Scalar& gradWSPH,
                      typename Dimension::Vector& gradWCSPH);

}
}

#ifndef __GCCXML__
#include "CSPHUtilitiesInline.hh"
#endif

#endif
