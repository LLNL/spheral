text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "CullenDehnenViscosity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class CullenDehnenViscosity< Dim< %(ndim)s > >;
  }
}
"""
