//---------------------------------Spheral++----------------------------------//
// CSPHUtilities
//
// Useful methods for using the CSPH formalism.
//
// Created by JMO, Fri Aug  8 16:16:33 PDT 2008
//----------------------------------------------------------------------------//
#include "Kernel/TableKernel.hh"

namespace Spheral {
namespace CSPHSpace {

//------------------------------------------------------------------------------
// Compute the corrected kernel value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
CSPHKernel(const KernelSpace::TableKernel<Dimension>& W,
           const typename Dimension::Vector& rij,
           const typename Dimension::Vector& etaj,
           const typename Dimension::Scalar& Hdetj,
           const typename Dimension::Scalar& Ai,
           const typename Dimension::Vector& Bi) {
  return Ai*(1.0 + Bi.dot(rij))*W(etaj.magnitude(), Hdetj);
}

//------------------------------------------------------------------------------
// Compute the corrected kernel value and gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
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
                      typename Dimension::Vector& gradWCSPH) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  const std::pair<Scalar, Scalar> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
  const Scalar Wj = WWj.first;
  WCSPH = Ai*(1.0 + Bi.dot(rij))*Wj;
  gradWSPH = WWj.second;
  const Vector gradWj = Hj*etaj.unitVector() * gradWSPH;
  gradWCSPH = Ai*(1.0 + Bi.dot(rij))*gradWj + Ai*(Bi + gradBi*rij)*Wj + gradAi*(1.0 + Bi.dot(rij))*Wj;
}

}
}
