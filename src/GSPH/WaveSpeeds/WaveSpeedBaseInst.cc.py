text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/WaveSpeedBase.cc"

namespace Spheral {
  template class WaveSpeedBase<Dim< %(ndim)s > >;
}
"""
