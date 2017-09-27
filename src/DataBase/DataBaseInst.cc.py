text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "DataBase.cc"

namespace Spheral {
  namespace DataBaseSpace {
    template class DataBase< Dim< %(ndim)s > >;
  }
}
"""
