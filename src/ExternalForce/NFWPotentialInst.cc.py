text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ExternalForce/NFWPotential.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class NFWPotential<Dim< %(ndim)s > >;
}
"""
