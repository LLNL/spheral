//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "findMatchingVertex.cc"

namespace Spheral {
  namespace MeshSpace {
    template unsigned findMatchingVertex<Dim<2>::Vector>(const Dim<2>::Vector& target,
                                                         const std::vector<Dim<2>::Vector>& verticesj);
    template unsigned findMatchingVertex<Dim<2>::Vector>(const Dim<2>::Vector& target,
                                                         const std::vector<Dim<2>::Vector>& verticesj,
                                                         const std::vector<unsigned>& indicesj);

    template unsigned findMatchingVertex<Dim<3>::Vector>(const Dim<3>::Vector& target,
                                                         const std::vector<Dim<3>::Vector>& verticesj);
    template unsigned findMatchingVertex<Dim<3>::Vector>(const Dim<3>::Vector& target,
                                                         const std::vector<Dim<3>::Vector>& verticesj,
                                                         const std::vector<unsigned>& indicesj);
  }
}

