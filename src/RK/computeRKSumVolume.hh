//---------------------------------Spheral++------------------------------------
// Compute the RK mass density summation.
//------------------------------------------------------------------------------
#ifndef __Spheral__computeRKSumVolume__
#define __Spheral__computeRKSumVolume__

#include "Field/FieldList.hh"
#include "Neighbor/ConnectivityMap.hh"
#include "Kernel/TableKernel.hh"
#include "Boundary/Boundary.hh"

#include <vector>

namespace Spheral {

template<typename Dimension>
void
computeRKSumVolume(const ConnectivityMap<Dimension>& connectivityMap,
                   const TableKernel<Dimension>& W,
                   const FieldList<Dimension, typename Dimension::Vector>& position,
                   const FieldList<Dimension, typename Dimension::Scalar>& mass,
                   const FieldList<Dimension, typename Dimension::SymTensor>& H,
                   FieldList<Dimension, typename Dimension::Scalar>& vol);

}

#endif
