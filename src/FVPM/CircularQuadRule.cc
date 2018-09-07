//---------------------------------Spheral++----------------------------------//
// CircularQuadRule
//
// Created by JNJ, Sun Jul 11 12:50:51 PDT 2010
//----------------------------------------------------------------------------//

#include "CircularQuadRule.hh"

namespace Spheral {

//-------------------------------------------------------------------
template <typename Dimension>
CircularQuadRule<Dimension>::
CircularQuadRule(const TableKernel<Dimension>& W):
  QuadRule<Dimension>(W)
{
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template <typename Dimension>
CircularQuadRule<Dimension>::
~CircularQuadRule()
{
}
//-------------------------------------------------------------------

//-------------------------------------------------------------------
template <typename Dimension>
void 
CircularQuadRule<Dimension>::
setDomain(const Vector& x1, 
          const SymTensor& H1,
          const Vector& x2,
          const SymTensor& H2)
{
  // We assume that the smoothing tensors are isotropic.
  double detH1 = H1.det();
  double r1 = this->mW.kernelExtent() / Dimension::rootnu(detH1);
  double detH2 = H2.det();
  double r2 = this->mW.kernelExtent() / Dimension::rootnu(detH2);
  setDomain(x1, r1, x2, r2);
}
//-------------------------------------------------------------------

}
