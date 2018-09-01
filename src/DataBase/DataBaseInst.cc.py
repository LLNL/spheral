text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase.cc"

namespace Spheral {
  template class DataBase< Dim< %(ndim)s > >;
}
"""
