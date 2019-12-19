text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ArtificialViscosity/CullenDehnenViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CullenDehnenViscosity< Dim< %(ndim)s > >;
}
"""
