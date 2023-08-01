text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/Policies/PairwisePlasticStrainPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class PairwisePlasticStrainPolicy<Dim< %(ndim)s > >;
}

"""
