text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "invertRankNTensor.cc"
#include "GeomFourthRankTensor.hh"

namespace Spheral {
  template GeomFourthRankTensor<%(ndim)s> invertRankNTensor(const GeomFourthRankTensor<%(ndim)s>&);
}
"""
