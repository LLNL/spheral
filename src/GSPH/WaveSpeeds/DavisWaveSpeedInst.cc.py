text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/DavisWaveSpeed.cc"

namespace Spheral {
  template class DavisWaveSpeed<Dim< %(ndim)s > >;
}
"""
