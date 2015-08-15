text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeHVolumes.cc"

namespace Spheral {
  template void computeHVolumes(const Dim<1>::Scalar kernelExtent,
                                const FieldList<Dim<1>, Dim<1>::SymTensor>&,
                                FieldList<Dim<1>, Dim<1>::Scalar>&);
}

"""
