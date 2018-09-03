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
  FieldList<Dimension, DataType>
  computeShepardsInterpolation(const FieldList<Dimension, DataType>& fieldList,
                               const ConnectivityMap<Dimension>& connectivityMap,
                               const TableKernel<Dimension>& W,
                               const FieldList<Dimension, typename Dimension::Vector>& position,
                               const FieldList<Dimension, typename Dimension::SymTensor>& H,
                               const FieldList<Dimension, typename Dimension::Scalar>& weight);
}

#endif
