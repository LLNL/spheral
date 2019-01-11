text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ArtificialConduction.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ArtificialConduction< Dim< %(ndim)s > >;
}
"""
