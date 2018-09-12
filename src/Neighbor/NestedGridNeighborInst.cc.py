text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Neighbor/NestedGridNeighbor.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {

//------------------------------------------------------------------------------
// Initialize static variables.
//------------------------------------------------------------------------------
template<> const double NestedGridNeighbor< Dim< %(ndim)s> >::ln2inverse = 1.0/log(2.0);
template<> const int NestedGridNeighbor< Dim< %(ndim)s> >::mEndOfLinkList = -1;
template<> const int NestedGridNeighbor< Dim< %(ndim)s > >::mGridNormalMagnitude = 1024;

template class NestedGridNeighbor< Dim< %(ndim)s > >;

}
"""
