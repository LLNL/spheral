text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "PredictorCorrector.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PredictorCorrector< Dim< %(ndim)s > >;
}
"""
