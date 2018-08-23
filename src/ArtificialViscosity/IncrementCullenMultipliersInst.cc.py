text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "ArtificialViscosity/IncrementCullenMultipliers.cc"

namespace Spheral {
  template class IncrementCullenMultipliers<Dim< %(ndim)s > >;
}
"""
