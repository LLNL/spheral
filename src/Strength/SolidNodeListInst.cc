//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "SolidNodeList.cc"

namespace Spheral {
  namespace SolidMaterial {
    template class SolidNodeList< Dim<1> >;
    template class SolidNodeList< Dim<2> >;
    template class SolidNodeList< Dim<3> >;
  }
}
