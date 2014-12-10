// //------------------------------------------------------------------------------
// // Static initializations.
// //------------------------------------------------------------------------------
// template<typename Dimension> typename Dimension::Vector Neighbor<Dimension>::mXmin = Dimension::Vector::zero;
// template<typename Dimension> typename Dimension::Vector Neighbor<Dimension>::mXmax = Dimension::Vector::zero;
#include "Neighbor.cc"

//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
namespace Spheral {
  namespace NeighborSpace {
    template class Neighbor< Dim<1> >;
    template class Neighbor< Dim<2> >;
    template class Neighbor< Dim<3> >;
  }
}
