//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ConnectivityMap.cc"
#include "Geometry/Dimension.hh"

template class Spheral::NeighborSpace::ConnectivityMap<Spheral::Dim<1> >;
template class Spheral::NeighborSpace::ConnectivityMap<Spheral::Dim<2> >;
template class Spheral::NeighborSpace::ConnectivityMap<Spheral::Dim<3> >;
