//---------------------------------Spheral++------------------------------------
// Compute the Shepard interpolation of a FieldList.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeShepardsInterpolation__
#define __Spheral__computeShepardsInterpolation__

#include <vector>

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"

namespace Spheral {
  template<typename Dimension, typename DataType>
  FieldSpace::FieldList<Dimension, DataType>
  computeShepardsInterpolation(const FieldSpace::FieldList<Dimension, DataType>& fieldList,
                               const NeighborSpace::ConnectivityMap<Dimension>& connectivityMap,
                               const KernelSpace::TableKernel<Dimension>& W,
                               const FieldSpace::FieldList<Dimension, typename Dimension::Vector>& position,
                               const FieldSpace::FieldList<Dimension, typename Dimension::SymTensor>& H,
                               const FieldSpace::FieldList<Dimension, typename Dimension::Scalar>& weight);
}

#endif
