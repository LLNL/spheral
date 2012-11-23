//---------------------------------Spheral++----------------------------------//
// HatKernel -- The B spline interpolation kernel.
//
// Created by JMO, Wed Dec 11 17:33:57 PST 2002
//----------------------------------------------------------------------------//

#include <math.h>

#include "Kernel.hh"
#include "HatKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class HatKernel<Dim<1> >;
    template class HatKernel<Dim<2> >;
    template class HatKernel<Dim<3> >;
  }
}
