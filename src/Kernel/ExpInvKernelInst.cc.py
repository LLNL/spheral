text = """
#include "Kernel/Kernel.hh"
#include "Kernel/ExpInvKernel.hh"

#include <math.h>

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  template class ExpInvKernel< Dim< %(ndim)s >  >;
}
"""
