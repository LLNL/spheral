text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CollinsStrength.cc"
#include "Geometry/Dimension.hh"

using namespace Spheral::Material;

namespace Spheral {
  namespace SolidMaterial {
    template class CollinsStrength<Dim< %(ndim)s > >;
  }
}
"""
