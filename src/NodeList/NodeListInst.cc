//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList.cc"

namespace Spheral {
  namespace NodeSpace {
    template class NodeList< Dim<1> >;
    template class NodeList< Dim<2> >;
    template class NodeList< Dim<3> >;
  }
}
