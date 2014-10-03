//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ArtificialConduction.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
    namespace PhysicsSpace {
        template class ArtificialConduction< Dim<1> >;
        template class ArtificialConduction< Dim<2> >;
        template class ArtificialConduction< Dim<3> >;
    }
}