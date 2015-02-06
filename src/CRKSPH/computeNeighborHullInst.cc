//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeNeighborHull.cc"

namespace Spheral {
  template 
  Dim<1>::FacetedVolume
  computeNeighborHull(const std::vector<std::vector<int> >& fullConnectivity,
                      const Dim<1>::Scalar etaCutoff,
                      const Dim<1>::Vector& ri,
                      const Dim<1>::SymTensor& Hi,
                      const FieldSpace::FieldList<Dim<1>, Dim<1>::Vector>& position);

  template 
  Dim<2>::FacetedVolume
  computeNeighborHull(const std::vector<std::vector<int> >& fullConnectivity,
                      const Dim<2>::Scalar etaCutoff,
                      const Dim<2>::Vector& ri,
                      const Dim<2>::SymTensor& Hi,
                      const FieldSpace::FieldList<Dim<2>, Dim<2>::Vector>& position);

  template 
  Dim<3>::FacetedVolume
  computeNeighborHull(const std::vector<std::vector<int> >& fullConnectivity,
                      const Dim<3>::Scalar etaCutoff,
                      const Dim<3>::Vector& ri,
                      const Dim<3>::SymTensor& Hi,
                      const FieldSpace::FieldList<Dim<3>, Dim<3>::Vector>& position);
}

