//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NestedGridNeighbor.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  namespace NeighborSpace {
    template class NestedGridNeighbor< Dim<1> >;
    template class NestedGridNeighbor< Dim<2> >;
    template class NestedGridNeighbor< Dim<3> >;

    // template int NestedGridNeighbor< Dim<1> >::gridLevel(const Dim<1>::Scalar&) const;
    // template int NestedGridNeighbor< Dim<2> >::gridLevel(const Dim<2>::Scalar&) const;
    // template int NestedGridNeighbor< Dim<3> >::gridLevel(const Dim<3>::Scalar&) const;

    // template int NestedGridNeighbor< Dim<1> >::gridLevel(const Dim<1>::SymTensor&) const;
    // template int NestedGridNeighbor< Dim<2> >::gridLevel(const Dim<2>::SymTensor&) const;
    // template int NestedGridNeighbor< Dim<3> >::gridLevel(const Dim<3>::SymTensor&) const;
  }
}
