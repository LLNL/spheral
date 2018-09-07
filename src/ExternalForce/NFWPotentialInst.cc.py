text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NFWPotential.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class NFWPotential<Dim< %(ndim)s > >;
}
"""
