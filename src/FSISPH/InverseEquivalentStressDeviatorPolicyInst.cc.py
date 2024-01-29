text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "FSISPH/InverseEquivalentStressDeviatorPolicy.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {
  template class InverseEquivalentStressDeviatorPolicy<Dim< %(ndim)s > >;
}

"""
