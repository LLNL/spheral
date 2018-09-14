text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/invertRankNTensor.cc"
#include "Geometry/GeomFourthRankTensor.hh"

namespace Spheral {
  template GeomFourthRankTensor<%(ndim)s> invertRankNTensor(const GeomFourthRankTensor<%(ndim)s>&);
}
"""
