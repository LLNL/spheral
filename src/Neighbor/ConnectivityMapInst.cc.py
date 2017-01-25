text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ConnectivityMap.cc"
#include "Geometry/Dimension.hh"

template class Spheral::NeighborSpace::ConnectivityMap<Spheral::Dim< %(ndim)s > >;
"""
