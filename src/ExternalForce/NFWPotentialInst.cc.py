text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NFWPotential.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace PhysicsSpace {
    using namespace Material;
    template class NFWPotential<Dim< %(ndim)s > >;
  }
}
"""
