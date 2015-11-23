text = """
//---------------------------------Spheral++----------------------------------//
// WendlandC2Kernel -- .
//----------------------------------------------------------------------------//

#include <math.h>

#include "Kernel.hh"
#include "WendlandC2Kernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class WendlandC2Kernel< Dim< %(ndim)s >  >;
  }
}
"""
