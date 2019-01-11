text = """
//---------------------------------Spheral++----------------------------------//
// QuarticSplineKernel -- A quartic spline, as described in
// Belytschko et al., Computational Methods in Applied Mathematics and Engineering
// 1996, 139, 3-47.
//
// Kernel extent: 2.0
//
// Created by JMO, Wed Jan  8 22:45:10 PST 2003
//----------------------------------------------------------------------------//
#include "Kernel/QuarticSplineKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class QuarticSplineKernel< Dim< %(ndim)s >  >;
}
"""
