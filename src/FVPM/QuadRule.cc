//---------------------------------Spheral++----------------------------------//
// QuadRule
//
// Created by JNJ, Sun Jul 11 12:50:51 PDT 2010
//----------------------------------------------------------------------------//

#include "QuadRule.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

namespace FVPMSpace {

//-------------------------------------------------------------------
template <typename Dimension>
QuadRule<Dimension>::
QuadRule(const KernelSpace::TableKernel<Dimension>& W):
  mW(W)
{
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
template <typename Dimension>
QuadRule<Dimension>::
~QuadRule()
{
}
//-------------------------------------------------------------------

}
}

