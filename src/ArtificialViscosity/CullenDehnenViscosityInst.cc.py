text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "CullenDehnenViscosity.cc"

namespace Spheral {
  namespace ArtificialViscositySpace {
    template class CullenDehnenViscosity< Dim< %(ndim)s > >;
  }
}
"""
