//---------------------------------Spheral++------------------------------------
// Compute the occupancy volume per point
//------------------------------------------------------------------------------
#ifndef __Spheral__computeOccupancyVolume__
#define __Spheral__computeOccupancyVolume__

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
void
computeOccupancyVolume(const ConnectivityMap<Dimension>& connectivityMap,
                       const TableKernel<Dimension>& W,
                       const FieldList<Dimension, typename Dimension::Vector>& position,
                       const FieldList<Dimension, typename Dimension::SymTensor>& H,
                       FieldList<Dimension, typename Dimension::Scalar>& vol);

}

#endif
