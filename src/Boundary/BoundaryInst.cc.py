text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary.cc"

namespace Spheral {
namespace BoundarySpace {
template class Boundary< Dim< %(ndim)s > >;
}
}

"""
