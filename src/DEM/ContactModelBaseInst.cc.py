text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "DEM/ContactModelBase.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class ContactModelBase< Dim< %(ndim)s > >;
}
"""
