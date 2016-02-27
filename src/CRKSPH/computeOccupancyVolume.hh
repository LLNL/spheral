//---------------------------------Spheral++------------------------------------
// Compute the occupancy volume per point
//------------------------------------------------------------------------------
#ifndef __Spheral__computeOccupancyVolume__
#define __Spheral__computeOccupancyVolume__

#include <vector>

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {
  namespace CRKSPHSpace {

    template<typename Dimension>
    void
    computeOccupancyVolume(const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                           const KernelSpace::TableKernel<Dimension>& W,
                           const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                           const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                           FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& vol);
  }
}

#endif
