text = """
//---------------------------------Spheral++----------------------------------//
// HatKernel -- The B spline interpolation kernel.
//
// Created by JMO, Wed Dec 11 17:33:57 PST 2002
//----------------------------------------------------------------------------//
#include "Kernel/Kernel.hh"
#include "Kernel/HatKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class HatKernel< Dim< %(ndim)s >  >;
}
"""
