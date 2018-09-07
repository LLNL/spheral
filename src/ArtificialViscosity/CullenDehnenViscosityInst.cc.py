text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CullenDehnenViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class CullenDehnenViscosity< Dim< %(ndim)s > >;
}
"""
