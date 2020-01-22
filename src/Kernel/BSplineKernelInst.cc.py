text = """
//---------------------------------Spheral++----------------------------------//
// BSplineKernel -- The B spline interpolation kernel.
//
// Created by JMO, Mon Nov 29 22:57:26 PST 1999
//----------------------------------------------------------------------------//
#include "Kernel/Kernel.hh"
#include "Kernel/BSplineKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class BSplineKernel<Dim< %(ndim)s > >;
}
"""
