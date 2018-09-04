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
#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Compute the corrected kernel value.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
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
             const typename Dimension::Scalar correctionMin,
             const typename Dimension::Scalar correctionMax) {
  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

  const Scalar Wij = 0.5*(W(etai.magnitude(), Hdeti) + W(etaj.magnitude(), Hdetj));
  if (correctionOrder == CRKOrder::ZerothOrder) {
    return std::max(correctionMin, 
                    std::min(correctionMax,
                             Ai))*Wij;
  } else if (correctionOrder == CRKOrder::LinearOrder) {
    return std::max(correctionMin,
                    std::min(correctionMax, 
                             Ai*(1.0 + Bi.dot(rij))))*Wij;
  } else {   //correctionOrder == QuadraticOrder
    return std::max(correctionMin, 
                    std::min(correctionMax, 
                             Ai*(1.0 + Bi.dot(rij) + innerDoubleProduct<Dimension>(Ci, rij.selfdyad()))))*Wij;
  }
}

//------------------------------------------------------------------------------
// Compute the corrected kernel value and gradient.
//------------------------------------------------------------------------------
template<typename Dimension>
inline
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
                        const typename Dimension::Scalar correctionMin,
                        const typename Dimension::Scalar correctionMax) {

  typedef typename Dimension::Scalar Scalar;
  typedef typename Dimension::Vector Vector;
  typedef typename Dimension::Tensor Tensor;

  // ij
  const std::pair<Scalar, Scalar> WWi = W.kernelAndGradValue(etai.magnitude(), Hdeti);
  const std::pair<Scalar, Scalar> WWj = W.kernelAndGradValue(etaj.magnitude(), Hdetj);
  const Scalar Wij = 0.5*(WWi.first + WWj.first); 
  const Vector gradWij = 0.5*(Hi*etai.unitVector() * WWi.second + Hj*etaj.unitVector() * WWj.second);
  gradWSPH = 0.5*(WWi.second + WWj.second);

  if (correctionOrder == CRKOrder::ZerothOrder) {
    const double correction0 = Ai;
    const double correction = std::max(correctionMin, std::min(correctionMax, correction0));
    WCRKSPH = correction*Wij;
    gradWCRKSPH = correction*gradWij;
    if (correction0 > correctionMin and correction0 < correctionMax) {
      gradWCRKSPH += gradAi*Wij;
    }

  } else if (correctionOrder == CRKOrder::LinearOrder) {
    const double correction0 = Ai*(1.0 + Bi.dot(rij));
    const double correction = std::max(correctionMin, std::min(correctionMax, correction0));
    WCRKSPH = correction*Wij;
    gradWCRKSPH = correction*gradWij;
    if (correction0 > correctionMin and correction0 < correctionMax) {
      gradWCRKSPH += Ai*Bi*Wij + gradAi*(1.0 + Bi.dot(rij))*Wij;
      for (size_t ii = 0; ii != Dimension::nDim; ++ii) {
        for (size_t jj = 0; jj != Dimension::nDim; ++jj) {
          gradWCRKSPH(ii) += Ai*Wij*gradBi(jj,ii)*rij(jj);
        }
      }
    }

  } else {  //correctionOrder == CRKOrder::QuadraticOrder
    const double correction0 = Ai*(1.0 + Bi.dot(rij) + innerDoubleProduct<Dimension>(Ci, rij.selfdyad()));
    const double correction = std::max(correctionMin, std::min(correctionMax, correction0));
    WCRKSPH = correction*Wij;
    gradWCRKSPH = correction*gradWij;
    if (correction0 > correctionMin and correction0 < correctionMax) {
      gradWCRKSPH += Ai*Bi*Wij;
      gradWCRKSPH += gradAi*(1.0 + Bi.dot(rij) + innerDoubleProduct<Dimension>(Ci, rij.selfdyad()))*Wij;
      gradWCRKSPH += Ai*(innerProduct<Dimension>(rij,gradBi))*Wij;
      gradWCRKSPH += Ai*(innerDoubleProduct<Dimension>(rij.selfdyad(),gradCi))*Wij;
      gradWCRKSPH += 2.0*Ai*(innerProduct<Dimension>(rij,Ci))*Wij;
    }
  }
}

}
