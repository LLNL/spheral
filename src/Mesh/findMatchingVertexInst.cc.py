text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "findMatchingVertexTempl.cc"

namespace Spheral {
  namespace MeshSpace {
    template unsigned findMatchingVertex<Dim< %(ndim)s >::Vector>(const Dim< %(ndim)s >::Vector& target,
                                                                  const std::vector<Dim< %(ndim)s >::Vector>& verticesj);
    template unsigned findMatchingVertex<Dim< %(ndim)s >::Vector>(const Dim< %(ndim)s >::Vector& target,
                                                                  const std::vector<Dim< %(ndim)s >::Vector>& verticesj,
                                                                  const std::vector<unsigned>& indicesj);
  }
}
"""
