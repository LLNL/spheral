text = """
//---------------------------------Spheral++----------------------------------//
// WendlandC6Kernel -- .
//
// Created by CDR, Nov 19 2014
//----------------------------------------------------------------------------//

#include <math.h>

#include "Kernel/Kernel.hh"
#include "Kernel/WendlandC6Kernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class WendlandC6Kernel< Dim< %(ndim)s >  >;
  }
}
"""
