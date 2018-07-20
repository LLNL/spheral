text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "StrengthModel.cc"
#include "Geometry/Dimension.hh"

using namespace Spheral::Material;

namespace Spheral {
  namespace SolidMaterial {
    template class StrengthModel<Dim< %(ndim)s > >;
  }
}
"""
