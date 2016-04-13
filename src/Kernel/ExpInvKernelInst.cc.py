text = """
#include <math.h>

#include "Kernel.hh"
#include "ExpInvKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class ExpInvKernel< Dim< %(ndim)s >  >;
  }
}
"""
