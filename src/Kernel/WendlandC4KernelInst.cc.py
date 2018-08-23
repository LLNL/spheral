text = """
//---------------------------------Spheral++----------------------------------//
// WendlandC4Kernel -- .
//
// Created by CDR, Nov 5 2014
//----------------------------------------------------------------------------//

#include <math.h>

#include "Kernel/Kernel.hh"
#include "Kernel/WendlandC4Kernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class WendlandC4Kernel< Dim< %(ndim)s >  >;
  }
}
"""
