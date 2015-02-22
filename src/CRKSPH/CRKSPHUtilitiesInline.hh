//---------------------------------Spheral++----------------------------------//
// CRKSPHUtilities
//
// Useful methods for using the CRKSPH formalism.
//
// Created by JMO, Fri Aug  8 16:16:33 PDT 2008
//----------------------------------------------------------------------------//
#include "Kernel/TableKernel.hh"

namespace Spheral {
namespace CRKSPHSpace {

//------------------------------------------------------------------------------
// Compute the corrected kernel value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
CRKSPHKernel(const KernelSpace::TableKernel<Dimension>& W,
           const typename Dimension::Vector& rij,
           const typename Dimension::Vector& etai,
           const typename Dimension::Scalar& Hdeti,
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
                      typename Dimension::Vector& gradWCRKSPH) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  const std::pair<Scalar, Scalar> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
  const Scalar Wj = WWj.first;
  WCRKSPH = Ai*(1.0 + Bi.dot(rij))*Wj;
  gradWSPH = WWj.second;
  const Vector gradWj = Hj*etaj.unitVector() * WWj.second;
  //gradWCRKSPH = Ai*(1.0 + Bi.dot(rij))*gradWj + Ai*(Bi + gradBi*rij)*Wj + gradAi*(1.0 + Bi.dot(rij))*Wj;
  gradWCRKSPH = Ai*(1.0 + Bi.dot(rij))*gradWj + Ai*Bi*Wj + gradAi*(1.0 + Bi.dot(rij))*Wj;
  for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
    for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
      gradWCRKSPH(ii) += Ai*Wj*gradBi(jj,ii)*rij(jj);
    }
  }
}

}
}
