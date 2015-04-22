#include "ArtificialConductionPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
    template class ArtificialConductionPolicy< Dim<1> >;
    template class ArtificialConductionPolicy< Dim<2> >;
    template class ArtificialConductionPolicy< Dim<3> >;
}