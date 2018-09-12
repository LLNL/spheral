text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Mesh/findMatchingVertexTempl.cc"

namespace Spheral {
  template unsigned findMatchingVertex<Dim< %(ndim)s >::Vector>(const Dim< %(ndim)s >::Vector& target,
                                                                const std::vector<Dim< %(ndim)s >::Vector>& verticesj);
  template unsigned findMatchingVertex<Dim< %(ndim)s >::Vector>(const Dim< %(ndim)s >::Vector& target,
                                                                const std::vector<Dim< %(ndim)s >::Vector>& verticesj,
                                                                const std::vector<unsigned>& indicesj);
}
"""
