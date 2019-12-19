text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Integrator/PredictorCorrector.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PredictorCorrector< Dim< %(ndim)s > >;
}
"""
