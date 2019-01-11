text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary.cc"

namespace Spheral {
template class Boundary< Dim< %(ndim)s > >;
}

"""
