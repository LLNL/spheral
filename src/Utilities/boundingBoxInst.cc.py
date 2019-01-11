text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Utilities/boundingBox.cc"

namespace Spheral {
  template void boundingBox(const vector<Dim< %(ndim)s >::Vector>& positions,
                            Dim< %(ndim)s >::Vector& xmin,
                            Dim< %(ndim)s >::Vector& xmax);

  template void boundingBox(const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& positions,
                            Dim< %(ndim)s >::Vector& xmin,
                            Dim< %(ndim)s >::Vector& xmax,
                            const bool useGhosts);
}
"""
