//---------------------------------Spheral++----------------------------------//
// CRKSPHUtilities
//
// Useful methods for using the CRKSPH formalism.
//
// Created by JMO, Fri Aug  8 16:16:33 PDT 2008
//----------------------------------------------------------------------------//
#include "Kernel/TableKernel.hh"
#include "Geometry/innerDoubleProduct.hh"
#include "Geometry/innerProduct.hh"
#include "CRKSPHCorrectionParams.hh"

namespace Spheral {
namespace CRKSPHSpace {

//------------------------------------------------------------------------------
// Compute the corrected kernel value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
typename Dimension::Scalar
CRKSPHKernel(const KernelSpace::TableKernel<Dimension>& W,
           const CRKOrder correctionOrder,
           const typename Dimension::Vector& rij,
           const typename Dimension::Vector& etai,
           const typename Dimension::Scalar& Hdeti,
           const typename Dimension::Vector& etaj,
           const typename Dimension::Scalar& Hdetj,
           const typename Dimension::Scalar& Ai,
           const typename Dimension::Vector& Bi,
           const typename Dimension::Tensor& Ci) {
  typedef typename Dimension::Tensor Tensor;
  if(correctionOrder == ZerothOrder){
      return Ai*W(etaj.magnitude(), Hdetj);
  }else if(correctionOrder == LinearOrder){
      return Ai*(1.0 + Bi.dot(rij))*W(etaj.magnitude(), Hdetj);
  }else {//correctionOrder == QuadraticOrder
    return Ai*(1.0 + Bi.dot(rij) + Geometry::innerDoubleProduct<Dimension>(Ci, rij.selfdyad()))*W(etaj.magnitude(), Hdetj);
  }
}

//------------------------------------------------------------------------------
// Compute the corrected kernel value and gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
void
CRKSPHKernelAndGradient(const KernelSpace::TableKernel<Dimension>& W,
                      const CRKOrder correctionOrder,
                      const typename Dimension::Vector& rij,
                      const typename Dimension::Vector& etai,
                      const typename Dimension::SymTensor& Hi,
                      const typename Dimension::Scalar& Hdeti,
                      const typename Dimension::Vector& etaj,
                      const typename Dimension::SymTensor& Hj,
                      const typename Dimension::Scalar& Hdetj,
                      const typename Dimension::Scalar& Ai,
                      const typename Dimension::Vector& Bi,
                      const typename Dimension::Tensor& Ci,
                      const typename Dimension::Vector& gradAi,
                      const typename Dimension::Tensor& gradBi,
                      const typename Dimension::ThirdRankTensor& gradCi,
                      typename Dimension::Scalar& WCRKSPH,
                      typename Dimension::Scalar& gradWSPH,
                      typename Dimension::Vector& gradWCRKSPH) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;
  const std::pair<Scalar, Scalar> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
  const Scalar Wj = WWj.first; 
  if(correctionOrder == ZerothOrder){
     WCRKSPH = Ai*Wj;
     gradWSPH = WWj.second;
     const Vector gradWj = Hj*etaj.unitVector() * WWj.second;
     gradWCRKSPH = Ai*gradWj + gradAi*Wj;
  }else if(correctionOrder == LinearOrder){
     WCRKSPH = Ai*(1.0 + Bi.dot(rij))*Wj;
     gradWSPH = WWj.second;
     const Vector gradWj = Hj*etaj.unitVector() * WWj.second;
     gradWCRKSPH = Ai*(1.0 + Bi.dot(rij))*gradWj + Ai*Bi*Wj + gradAi*(1.0 + Bi.dot(rij))*Wj;
     for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
       for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
         gradWCRKSPH(ii) += Ai*Wj*gradBi(jj,ii)*rij(jj);
       }
     }
  }else {//correctionOrder == QuadraticOrder
     WCRKSPH = Ai*(1.0 + Bi.dot(rij) + Geometry::innerDoubleProduct<Dimension>(Ci, rij.selfdyad()))*Wj;
     gradWSPH = WWj.second;
     const Vector gradWj = Hj*etaj.unitVector() * WWj.second;
     gradWCRKSPH = Ai*(1.0 + Bi.dot(rij) + Geometry::innerDoubleProduct<Dimension>(Ci, rij.selfdyad()))*gradWj + Ai*Bi*Wj;
     gradWCRKSPH += gradAi*(1.0 + Bi.dot(rij) + Geometry::innerDoubleProduct<Dimension>(Ci, rij.selfdyad()))*Wj;
     gradWCRKSPH += Ai*(Geometry::innerProduct<Dimension>(rij,gradBi))*Wj;
     gradWCRKSPH += Ai*(Geometry::innerDoubleProduct<Dimension>(rij.selfdyad(),gradCi))*Wj;
     gradWCRKSPH += 2.0*Ai*(Geometry::innerProduct<Dimension>(rij,Ci))*Wj;
  }
}

}
}
