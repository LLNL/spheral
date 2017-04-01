text = """
//------------------------------------------------------------------------------
// Explict instantiation.
//------------------------------------------------------------------------------
#include "ArtificialConduction.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
    namespace PhysicsSpace {
        template class ArtificialConduction< Dim< %(ndim)s > >;
    }
}
"""
