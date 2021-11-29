text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/EinfeldtWaveSpeed.cc"

namespace Spheral {
  template class EinfeldtWaveSpeed<Dim< %(ndim)s > >;
}
"""
