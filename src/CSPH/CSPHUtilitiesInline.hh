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
           const typename Dimension::Vector& etai,
           const typename Dimension::Scalar& Hdeti,
           const typename Dimension::Vector& etaj,
           const typename Dimension::Scalar& Hdetj,
           const typename Dimension::Scalar& Ai,
           const typename Dimension::Vector& Bi) {
  return Ai*(1.0 + Bi.dot(rij))*0.5*(W(etaj.magnitude(), Hdetj) + 
                                     W(etai.magnitude(), Hdeti));
}

//------------------------------------------------------------------------------
// Compute the corrected kernel value and gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
CSPHKernelAndGradient(const KernelSpace::TableKernel<Dimension>& W,
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
                      typename Dimension::Scalar& WCSPH,
                      typename Dimension::Scalar& gradWSPH,
                      typename Dimension::Vector& gradWCSPH) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  const std::pair<Scalar, Scalar> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
  const std::pair<Scalar, Scalar> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
  const Scalar Wj = WWj.first;
  const Scalar Wi = WWi.first;
  const Scalar Wij = 0.5*(Wj + Wi);
  WCSPH = Ai*(1.0 + Bi.dot(rij))*Wij;
  gradWSPH = 0.5*(WWj.second + WWi.second);
  const Vector gradWij = 0.5*(Hj*etaj.unitVector() * WWj.second -
                              Hi*etai.unitVector() * WWi.second);
  //gradWCSPH = Ai*(1.0 + Bi.dot(rij))*gradWj + Ai*(Bi + gradBi*rij)*Wj + gradAi*(1.0 + Bi.dot(rij))*Wj;
  gradWCSPH = Ai*(1.0 + Bi.dot(rij))*gradWij + Ai*Bi*Wij + gradAi*(1.0 + Bi.dot(rij))*Wij;
  for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
    for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
      gradWCSPH(ii) += Ai*Wij*gradBi(jj,ii)*rij(jj);
    }
  }
}

}
}
