text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "GSPH/computeMFMDensity.cc"
#include "Geometry/Dimension.hh"

namespace Spheral {


  template void computeMFMDensity(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                  const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&,
                                        FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Scalar>&);
}
"""
