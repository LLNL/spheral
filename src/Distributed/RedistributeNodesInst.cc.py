text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Distributed/RedistributeNodes.cc"

namespace Spheral {
  template class RedistributeNodes< Dim< %(ndim)s > >;
}

"""
