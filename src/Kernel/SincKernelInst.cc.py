text = """
//---------------------------------Spheral++----------------------------------//
// SincKernel -- The sinc interpolation kernel.
//
// Created by JMO, Wed Dec  1 14:38:51 PST 1999
//----------------------------------------------------------------------------//
#include "Kernel/SincKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class SincKernel< Dim< %(ndim)s >  >;
}
"""
