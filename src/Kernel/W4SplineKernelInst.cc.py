text = """
//---------------------------------Spheral++----------------------------------//
// W4SplineKernel -- The B spline interpolation kernel.
//
// Created by JMO, Mon Nov 29 22:57:26 PST 1999
//----------------------------------------------------------------------------//

#include <math.h>

#include "Kernel.hh"
#include "W4SplineKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class W4SplineKernel< %(ndim)s  >;
  }
}
"""
