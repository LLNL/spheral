text = """
#include <math.h>

#include "Kernel/Kernel.hh"
#include "Kernel/ExpInvKernel.hh"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace KernelSpace {
    template class ExpInvKernel< Dim< %(ndim)s >  >;
  }
}
"""
