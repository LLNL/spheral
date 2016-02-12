text = """
//---------------------------------Spheral++----------------------------------//
// BSplineKernel -- The B spline interpolation kernel.
//
// Created by JMO, Mon Nov 29 22:57:26 PST 1999
//----------------------------------------------------------------------------//

#include <math.h>

#include "Kernel.hh"
#include "BSplineKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class BSplineKernel<Dim< %(ndim)s > >;
  }
}
"""
