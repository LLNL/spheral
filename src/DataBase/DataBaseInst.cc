//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase.cc"

namespace Spheral {
  namespace DataBaseSpace {
    template class DataBase< Dim<1> >;
    template class DataBase< Dim<2> >;
    template class DataBase< Dim<3> >;
  }
}
