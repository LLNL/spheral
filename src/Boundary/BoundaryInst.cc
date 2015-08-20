//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "Boundary.cc"

namespace Spheral {
namespace BoundarySpace {
template class Boundary< Dim<1> >;
template class Boundary< Dim<2> >;
template class Boundary< Dim<3> >;

}
}

