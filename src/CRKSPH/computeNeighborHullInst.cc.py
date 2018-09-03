text = """
//------------------------------------------------------------------------------
// Explicit instantiation.
//------------------------------------------------------------------------------
#include "Geometry/Dimension.hh"
#include "computeNeighborHull.cc"

namespace Spheral {
  template 
  Dim< %(ndim)s >::FacetedVolume
  computeNeighborHull(const std::vector<std::vector<int> >& fullConnectivity,
                      const Dim< %(ndim)s >::Scalar etaCutoff,
                      const Dim< %(ndim)s >::Vector& ri,
                      const Dim< %(ndim)s >::SymTensor& Hi,
                      const FieldList<Dim< %(ndim)s >, Dim< %(ndim)s >::Vector>& position);
}

"""
