//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "Mesh/findMatchingVertexTempl.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template unsigned findMatchingVertex<Dim<1>::Vector>(const Dim<1>::Vector& target,
                                                                const std::vector<Dim<1>::Vector>& verticesj);
  template unsigned findMatchingVertex<Dim<1>::Vector>(const Dim<1>::Vector& target,
                                                                const std::vector<Dim<1>::Vector>& verticesj,
                                                                const std::vector<unsigned>& indicesj);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template unsigned findMatchingVertex<Dim<2>::Vector>(const Dim<2>::Vector& target,
                                                                const std::vector<Dim<2>::Vector>& verticesj);
  template unsigned findMatchingVertex<Dim<2>::Vector>(const Dim<2>::Vector& target,
                                                                const std::vector<Dim<2>::Vector>& verticesj,
                                                                const std::vector<unsigned>& indicesj);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template unsigned findMatchingVertex<Dim<3>::Vector>(const Dim<3>::Vector& target,
                                                                const std::vector<Dim<3>::Vector>& verticesj);
  template unsigned findMatchingVertex<Dim<3>::Vector>(const Dim<3>::Vector& target,
                                                                const std::vector<Dim<3>::Vector>& verticesj,
                                                                const std::vector<unsigned>& indicesj);
#endif
}