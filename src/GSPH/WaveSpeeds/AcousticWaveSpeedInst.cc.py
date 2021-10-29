text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "GSPH/WaveSpeeds/AcousticWaveSpeed.cc"

namespace Spheral {
  template class AcousticWaveSpeed<Dim< %(ndim)s > >;
}
"""
