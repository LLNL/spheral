text = """
//---------------------------------Spheral++----------------------------------//
// WendlandC2Kernel -- .
//----------------------------------------------------------------------------//
#include "Kernel/Kernel.hh"
#include "Kernel/WendlandC2Kernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class WendlandC2Kernel< Dim< %(ndim)s >  >;
}
"""
