//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "ASPHSmoothingScale.cc"

namespace Spheral {
  namespace NodeSpace {
    template class ASPHSmoothingScale<Dim<1> >;
    template class ASPHSmoothingScale<Dim<2> >;
    template class ASPHSmoothingScale<Dim<3> >;
  }
}
