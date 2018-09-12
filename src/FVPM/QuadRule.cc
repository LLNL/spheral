//---------------------------------Spheral++----------------------------------//
// QuadRule
//
// Created by JNJ, Sun Jul 11 12:50:51 PDT 2010
//----------------------------------------------------------------------------//

#include "QuadRule.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {

//-------------------------------------------------------------------
template <typename Dimension>
QuadRule<Dimension>::
QuadRule(const TableKernel<Dimension>& W):
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
