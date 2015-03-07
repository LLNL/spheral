//---------------------------------Spheral++------------------------------------
// Compute the hull for a given points neighbor set.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeNeighborHull__
#define __Spheral__computeNeighborHull__

#include <vector>

namespace Spheral {

  // Forward declarations.
  namespace FieldSpace {
    template<typename Dimension, typename DataType> class FieldList;
  }

  template<typename Dimension>
  typename Dimension::FacetedVolume
  computeNeighborHull(const std::vector<std::vector<int> >& fullConnectivity,
                      const typename Dimension::Scalar etaCutoff,
                      const typename Dimension::Vector& ri,
                      const typename Dimension::SymTensor& Hi,
                      const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position);
}

#endif
