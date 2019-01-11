text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "NodeList/generateVoidNodes.cc"

namespace Spheral {
  template void generateVoidNodes<Dim< %(ndim)s > >(const vector<Dim< %(ndim)s >::Vector>& generators,
                                           const vector<Dim< %(ndim)s >::SymTensor>& Hs,
                                           const Mesh<Dim< %(ndim)s > >& mesh,
                                           const Dim< %(ndim)s >::Vector& xmin,
                                           const Dim< %(ndim)s >::Vector& xmax,
                                           const unsigned numInternal,
                                           const double nPerh,
                                           const double threshold,
                                           NodeList<Dim< %(ndim)s > >& voidNodes);
}
"""
