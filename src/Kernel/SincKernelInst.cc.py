text = """
//---------------------------------Spheral++----------------------------------//
// SincKernel -- The sinc interpolation kernel.
//
// Created by JMO, Wed Dec  1 14:38:51 PST 1999
//----------------------------------------------------------------------------//

#include <math.h>

#include "SincKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class SincKernel< Dim< %(ndim)s >  >;
  }
}
"""
