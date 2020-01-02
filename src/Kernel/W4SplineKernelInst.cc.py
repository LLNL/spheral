text = """
//---------------------------------Spheral++----------------------------------//
// W4SplineKernel -- The B spline interpolation kernel.
//
// Created by JMO, Mon Nov 29 22:57:26 PST 1999
//----------------------------------------------------------------------------//
#include "Kernel/Kernel.hh"
#include "Kernel/W4SplineKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class W4SplineKernel< Dim< %(ndim)s >  >;
}
"""
