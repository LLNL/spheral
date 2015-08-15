text = """
//---------------------------------Spheral++----------------------------------//
// QuinticSplineKernel -- A quintic spline, as described in
// Bonet & Kulasegaruam 2002, Appl. Math. Comput., 126, 133-155.
//
// Kernel extent: 2.0
//
// Created by JMO, Wed Jul  9 16:24:25 PDT 2008
//----------------------------------------------------------------------------//
#include "QuinticSplineKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class QuinticSplineKernel< %(ndim)s  >;
  }
}
"""
