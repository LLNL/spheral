text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "IncrementCullenMultipliers.cc"

namespace Spheral {
  template class IncrementCullenMultipliers<Dim< %(ndim)s > >;
}
"""
