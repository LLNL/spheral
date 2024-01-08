text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Porosity/PorositySolidMassDensityPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PorositySolidMassDensityPolicy<Dim< %(ndim)s > >;
}

"""
