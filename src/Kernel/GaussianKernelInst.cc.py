text = """
//---------------------------------Spheral++----------------------------------//
// GaussianKernel -- The gaussian interpolation kernel.
//
// Created by JMO, Wed Dec  1 14:38:51 PST 1999
//----------------------------------------------------------------------------//
#include "Kernel/Kernel.hh"
#include "Kernel/GaussianKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class GaussianKernel< Dim< %(ndim)s >  >;
}
"""
