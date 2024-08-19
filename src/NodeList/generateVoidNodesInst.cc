//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------

#include "config.hh"
#include "NodeList/generateVoidNodes.cc"

namespace Spheral {
#if defined(SPHERAL_ENABLE_1D)
  template void generateVoidNodes<Dim<1> >(const vector<Dim<1>::Vector>& generators,
                                           const vector<Dim<1>::SymTensor>& Hs,
                                           const Mesh<Dim<1> >& mesh,
                                           const Dim<1>::Vector& xmin,
                                           const Dim<1>::Vector& xmax,
                                           const unsigned numInternal,
                                           const double nPerh,
                                           const double threshold,
                                           NodeList<Dim<1> >& voidNodes);
#endif

#if defined(SPHERAL_ENABLE_2D)
  template void generateVoidNodes<Dim<2> >(const vector<Dim<2>::Vector>& generators,
                                           const vector<Dim<2>::SymTensor>& Hs,
                                           const Mesh<Dim<2> >& mesh,
                                           const Dim<2>::Vector& xmin,
                                           const Dim<2>::Vector& xmax,
                                           const unsigned numInternal,
                                           const double nPerh,
                                           const double threshold,
                                           NodeList<Dim<2> >& voidNodes);
#endif

#if defined(SPHERAL_ENABLE_3D)
  template void generateVoidNodes<Dim<3> >(const vector<Dim<3>::Vector>& generators,
                                           const vector<Dim<3>::SymTensor>& Hs,
                                           const Mesh<Dim<3> >& mesh,
                                           const Dim<3>::Vector& xmin,
                                           const Dim<3>::Vector& xmax,
                                           const unsigned numInternal,
                                           const double nPerh,
                                           const double threshold,
                                           NodeList<Dim<3> >& voidNodes);
#endif
}