//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "TableKernel.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace KernelSpace {
    template class TableKernel< Dim<1> >;
    template class TableKernel< Dim<2> >;
    template class TableKernel< Dim<3> >;
  }
}
