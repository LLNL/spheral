text = """
#include "ArtificialConductionPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
    template class ArtificialConductionPolicy< Dim< %(ndim)s > >;
}
"""
