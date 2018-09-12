text = """
#include "Kernel.hh"
#include "ExpInvKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class ExpInvKernel< Dim< %(ndim)s >  >;
}
"""
