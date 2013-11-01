//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FluidNodeList.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NodeSpace {
    template class FluidNodeList< Dim<1> >;
    template class FluidNodeList< Dim<2> >;
    template class FluidNodeList< Dim<3> >;
  }
}
